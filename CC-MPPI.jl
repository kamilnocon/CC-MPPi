using Plots
using LinearAlgebra
using Distributions
include("parameters.jl")
include("VehicleModels.jl")
include("utils.jl")

# --------[Define Variables]-----------
x_dim = 4
u_dim = 2
N = 15
dt = 0.01

# MPPI Hyperparams
lambda = 1
M = 1024
mu = zeros(Float64, u_dim)

sigma_control = zeros(Float64, u_dim, u_dim)
sigma_control[1, 1] = 0.49 * 5
sigma_control[2, 2] = 0.12 * 5
sigma_xf = Matrix{Float64}(I, x_dim, x_dim)*0.001

# En = zeros(Float64, (x_dim, (N+1)*x_dim))
# En[:,N*x_dim+1:(N+1)*x_dim] = Matrix{Float64}(I, x_dim, x_dim)
En = zeros(Float64, (x_dim, (N)*x_dim))
En[:,(N-1)*x_dim+1:(N)*x_dim] = Matrix{Float64}(I, x_dim, x_dim)
Q = Matrix{Float64}(I, x_dim, x_dim)
R = Matrix{Float64}(I, u_dim, u_dim)*0.1
Q_bar = zeros(Float64, ((N)*x_dim, (N)*x_dim))
R_bar = zeros(Float64, (N*u_dim, N*u_dim))
sigma_control_bar = zeros(Float64, (N*u_dim, N*u_dim))
for k in 1:N
    Q_bar[(k-1)*x_dim + 1:k*x_dim , (k-1)*x_dim + 1:k*x_dim] = Q
    R_bar[(k-1)*u_dim + 1:k*u_dim , (k-1)*u_dim + 1:k*u_dim] = R
    sigma_control_bar[(k-1)*u_dim + 1:k*u_dim , (k-1)*u_dim + 1:k*u_dim] = sigma_control
end
Q_bar[(N-1)*x_dim + 1:N*x_dim , (N-1)*x_dim + 1:N*x_dim] = Q

function main()
    # set inital reference trajectory
    X_ref = zeros(Float64, N+1, x_dim)
    X_ref[1,:] = [-10.0, 0.0, 0.0, 5.0]
    U_ref = zeros(Float64, N, u_dim)

    X_log = zeros(Float64, 1, x_dim)
    X_log[1, :] = X_ref[1,:]

    # Keep track of all control squences and their costs for MPPI 
    Sm_list = zeros(Float64, M)
    Um_list = zeros(Float64, M, N, u_dim)

    # Set goals and obstacles
    goal = [10.0, 0.0]
    obs_info = [-2.0, 0.0, 1]
    obs_info2 = [2.0, 0.0, 1]

    # while task not complete
    i = 0
    # while abs(X_ref[1,1]-goal[1]) > 1
    anim = @animate for i in 1:50
        # roll out dynamics
        for k in 1:N
            X_ref[k+1,:] = BicycleModelDynamics(X_ref[k,:], U_ref[k,:])
        end
        # println(X_ref[1,:])
        print(i)
        i+=1

        A, B, K = CovarianceControl(X_ref, U_ref)

        for m in 1:M
            Sm = 0
            eps = rand(MvNormal(mu, sigma_control), N)
            Um = copy(U_ref)
            yk = zeros(Float64, x_dim)
            xk = X_ref[1,:]
            for k in 1:N
                Kk = K[(k-1)*u_dim + 1:k*u_dim , (k-1)*x_dim + 1:k*x_dim]
                if m < (1-0.2)*M
                    Um[k, :] += Kk * yk
                else
                    Um[k, :] = zeros(Float64, u_dim)
                end
                Um[k, :] += eps[:, k]
                xk = BicycleModelDynamics(xk, Um[k,:]) #update state
                yk = A[(k-1)*x_dim + 1:k*x_dim , (k-1)*x_dim + 1:k*x_dim] * yk + B[(k-1)*x_dim + 1:k*x_dim , (k-1)*u_dim + 1:k*u_dim] * eps[:, k]
                Sm = Sm + cost(xk, Um[k, :], obs_info, obs_info2, eps, R)
            end
            Sm = Sm + terminal_cost(xk, X_ref[1, 1])
            Sm_list[m] = Sm
            Um_list[m, :, :] = Um
        end

        #Calculate Optimal Control
        V = optimal_control(Sm_list, Um_list, lambda, M)
        #Execute Command
        X_ref[1,:] = BicycleModelDynamics(X_ref[1,:], V[1,:])'
        X_log = vcat(X_log, X_ref[1,:]')
        U_ref[1:N-1,:] = V[2:N,:]
        U_ref[N, :] = V[N, :]

        # plot test
        p = plot(size = [4000, 2000])
        for m in 1:M
            X_m = X_ref
            for k in 1:N
                X_m[k+1,:] = BicycleModelDynamics(X_m[k,:], Um_list[m,k,:])'
            end
            plot!(X_m[:,1], X_m[:,2], legend = false, ylims=(-5,5), xlims=(-10,10))
            phi = range(0, stop = 2*pi, length = 100)
            x_obs = obs_info[1] .+ cos.(phi)
            y_obs = obs_info[2] .+ sin.(phi)
            plot!(x_obs, y_obs, line = (7, :black))
            x_obs2 = obs_info2[1] .+ cos.(phi)
            y_obs2 = obs_info2[2] .+ sin.(phi)
            plot!(x_obs2, y_obs2, line = (7, :black))
            plot!(X_log[:,1], X_log[:,2], line = (7, :green))
        end
    end
    gif(anim, "test.gif", fps = 10)
end
main()