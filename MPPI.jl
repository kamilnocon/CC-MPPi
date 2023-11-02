using Plots
using LinearAlgebra
using Distributions
using Statistics
using DelimitedFiles
#include("parameters.jl")
include("VehicleModels.jl")
include("utils.jl")

# --------[Define Variables]-----------
x_dim = 4
u_dim = 2
N = 10
dt = 0.1
nu = 1000

# MPPI Hyperparams
lambda = 1500
M = 1024
mu = zeros(Float64, u_dim)

sigma_control = zeros(Float64, u_dim, u_dim)
sigma_control[1, 1] = 0.49
sigma_control[2, 2] = 0.12
sigma_xf = Matrix{Float64}(0.01I, x_dim, x_dim)
#=
sigma_xf[1, 1] = 0.1
sigma_xf[2, 2] = 0.1
sigma_xf[3, 3] = 0.01
sigma_xf[4, 4] = 0.01
=#

En = zeros(Float64, (x_dim, (N+1)*x_dim))
En[:,N*x_dim+1:(N+1)*x_dim] = Matrix{Float64}(I, x_dim, x_dim)
#Q = Matrix{Float64}(I, x_dim, x_dim)
Q = zeros(Float64, x_dim, x_dim)
R = Matrix{Float64}(I, u_dim, u_dim)*0.001
Q_bar = zeros(Float64, ((N+1)*x_dim, (N+1)*x_dim))
R_bar = zeros(Float64, (N*u_dim, N*u_dim))
sigma_control_bar = zeros(Float64, (N*u_dim, N*u_dim))
for k in 1:N
    Q_bar[(k-1)*x_dim + 1:k*x_dim , (k-1)*x_dim + 1:k*x_dim] = Q
    R_bar[(k-1)*u_dim + 1:k*u_dim , (k-1)*u_dim + 1:k*u_dim] = R
    sigma_control_bar[(k-1)*u_dim + 1:k*u_dim , (k-1)*u_dim + 1:k*u_dim] = sigma_control
end
Q_bar[N*x_dim + 1:(N+1)*x_dim , N*x_dim + 1:(N+1)*x_dim] = Q

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
    obs_info = [-7.5, 0.0, 0.75]
    obs_info2 = [4.0, 0.0, 0.75]

    # while task not complete
    i = 0
    # while abs(X_ref[1,1]-goal[1]) > 1
    anim = @animate for i in 1:30
    # anim = @animate while abs(X_ref[1,1]-goal[1]) > 1
        # roll out dynamics
        #=
        for k in 1:N
            X_ref[k+1,:] = BicycleModelDynamics(X_ref[k,:], U_ref[k,:])
        end
        =#
        print(i)
        i+=1

        #A, B, K = CovarianceControl(X_ref, U_ref)

        for m in 1:M
            Sm = 0
            eps = rand(MvNormal(mu, sigma_control), N)
            Um = copy(U_ref)
            yk = zeros(Float64, x_dim)
            xk = X_ref[1,:]
            for k in 1:N
                #=
                Kk = K[(k-1)*u_dim + 1:k*u_dim , (k)*x_dim + 1:(k+1)*x_dim]
                if m < (1-0.2)*M
                    Um[k, :] += Kk * yk
                else
                    Um[k, :] = zeros(Float64, u_dim)
                end
                =#
                Um[k, :] += eps[:, k]
                xk = BicycleModelDynamics(xk, Um[k,:]) #update state
                #yk = A[k,:,:] * yk + B[k,:,:] * eps[:, k]
                Sm += cost(xk, Um[k, :], obs_info, obs_info2, eps, R)
            end
            Sm += terminal_cost(xk, X_ref[1, 1])
            Sm_list[m] = Sm
            Um_list[m, :, :] = Um
        end

        X_best = copy(X_ref)

        #Calculate Optimal Control
        V = optimal_control(Sm_list, Um_list, lambda, M)
        #Execute Command
        #X_ref[1,:] = BicycleModelDynamics(X_ref[1,:], V[1,:])'
        X_ref[1,:] = ExecuteCommands(X_ref[1,:], V[1,:], dt)'
        X_log = vcat(X_log, X_ref[1,:]')
        U_ref[1:N-1,:] = V[2:N,:]
        U_ref[N, :] = V[N, :]

        # plot test
        p = plot(size = [4000, 2000])
        X_f_m = zeros(Float64, M, x_dim)
        for m in 1:M
            X_m = copy(X_ref)
            for k in 1:N
                X_m[k+1,:] = BicycleModelDynamics(X_m[k,:], Um_list[m,k,:])'
            end
            X_f_m[m,:] = X_m[N,:]
            plot!(X_m[:,1], X_m[:,2], line = (1, :grey), legend = false, ylims=(-5,5), xlims=(-10,10))
        end
        for k in 1:N
            X_best[k+1,:] = BicycleModelDynamics(X_best[k,:], V[k,:])'
        end
        plot!(X_best[:,1], X_best[:,2], line = (5, :green))
        phi = range(0, stop = 2*pi, length = 100)
        x_obs = obs_info[1] .+ obs_info[3]*cos.(phi)
        y_obs = obs_info[2] .+ obs_info[3]*sin.(phi)
        plot!(x_obs, y_obs, line = (7, :black))
        x_obs2 = obs_info2[1] .+ obs_info[3]*cos.(phi)
        y_obs2 = obs_info2[2] .+ obs_info[3]*sin.(phi)
        plot!(x_obs2, y_obs2, line = (7, :black))
        plot!(X_log[:,1], X_log[:,2], line = (7, :red))
        # println(round.(Statistics.cov(X_f_m, dims=1), digits=4))
    end
    gif(anim, "test.gif", fps = 10)
end
main()