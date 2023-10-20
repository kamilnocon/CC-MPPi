using Distributions
using JuMP
using Plots
using NPZ
using DelimitedFiles
using LinearAlgebra

function cost(x0, u0, obs_info, obs_info2, eps, R)
    state_cost = 0
    if ((x0[1] < -30.0) | (x0[1] > 30.0)) & ((x0[2] < -5.0) | (x0[2] > 5.0))
        state_cost = state_cost + 2000.0
    end
    if sqrt((x0[1]-obs_info[1])^2+(x0[2]-obs_info[2])^2) < obs_info[3]
        state_cost = state_cost + 10.0*712.5
    end
    if sqrt((x0[1]-obs_info2[1])^2+(x0[2]-obs_info2[2])^2) < obs_info2[3]
        state_cost = state_cost + 10.0*712.5
    end
    sampling_cost = (transpose(u0)*R*eps)[1] + (1/2)*(transpose(u0)*R*u0)[1]
    rollout_cost = state_cost + sampling_cost
    return rollout_cost
end

function terminal_cost(x0, x)
    s = abs(x0[1]-x)
    e = abs(0.0-x0[2])
    term_cost = 3.3*(1-s)+500*e^2
    return term_cost
end

function optimal_control(Sm_list, U0_list, lambda, M)
    b = minimum(Sm_list)
    V = zeros(Float64, length(U0_list[1, :]))
    sum = 0
    for m in 1:M
        wm = exp((-1/lambda)*(Sm_list[m]-b))
        V = V + wm.*U0_list[m, :]
        sum = sum + wm
    end
    V = V ./ sum
    return V
end

function dynamics_test(x0, u0, dt, lf, lr)
    beta = atan((lr/((lf+lr)))*tan(u0[2]))
    dx = zeros(Float64, 4)
    dx[1] = x0[4]*cos(beta+x0[3])
    dx[2] = x0[4]*sin(beta+x0[3])
    dx[3] = ((x0[4])*sin(beta))/(lr)
    dx[4] = u0[1]
    x0 = x0 + dt*dx
end



function main()
    N = 15
    dt = 0.01
    lambda = 1
    M = 1024
    numStates = 4
    numControls = 2
    R = Matrix(1.0I, numControls, numControls)

    X0 = zeros(Float64, N+1, numStates)
    U0 = zeros(Float64, N, numControls)
    #U0[1, :] = [0.0, 0.0]
    X = zeros(Float64, numStates)
    U = zeros(Float64, numControls)
    lf = 1
    lr = 1
    
    alpha = 0.2

    eps_cov = zeros(Float64, numControls, numControls)
    eps_cov[1, 1] = 0.0049
    eps_cov[2, 2] = 0.0012
    mu = zeros(Float64, numControls)

    Sm_list = zeros(Float64, M)
    U0_list = zeros(Float64, M, numControls*N)

    goal = [10.0, 0.0]
    obs_info = [-2.0, 0.0, 1]
    obs_info2 = [2.0, 0.0, 1]

    #X0[1, :] = [-25.0, 0.0, 0.0, 0.0]
    #eps = rand(MvNormal(mu, eps_cov), N)
    X = [-10.0, 0.0, 0.0, 5.0]
    X_record = zeros(Float64, 1, numStates)
    X_record[1, :] = X
    
    i = 0
    while abs(X[1]-goal[1]) > 1
        i+=1
        for m in 1:M
            # Um = U0
            Sm = 0
            #Y0 = Y0
            X0[1, :] = X
            eps = rand(MvNormal(mu, eps_cov), N)
            for k in 1:N
                U0[k, :] = U0[k, :] + eps[:, k] 
                for j in 1:numControls
                    U0_list[m, numControls*k-numControls+j] = U0[k, j]
                end
                X0[k+1, :] = dynamics_test(X0[k, :], U0[k, :], dt, lf, lr) #update state
                #Y0[k+1, :] = Y0 #update feedback
                Sm = Sm + cost(X0[k, :], U0[k, :], obs_info, obs_info2, eps, R)
            end
            Sm = Sm + terminal_cost(X0[N, :], X0[1, 1]) #add terminal cost here
            Sm_list[m] = Sm
        end
        #Calculate Optimal Control
        V = optimal_control(Sm_list, U0_list, lambda, M)
        U[1] = V[1]
        U[2] = V[2]
        #Execute Command
        X = dynamics_test(X, U, dt, lf, lr)
        #println(X)
        #println([X[1], X[2], maximum(Sm_list)])
        X_record = vcat(X_record, X')
        #Initialize W
        for i in 1:N-1
            U0[i, :] = getindex(V, numControls*(i+1)-numControls+1:numControls*(i+1)-numControls+2)
        end
        U0[N, :] = U0[N-1, :]
    end
    println(i)
    #println(X_record)
    plot!(X_record[:, 1], X_record[:, 2], tickfontsize=30, line = (7, :blue), xlabel = "X (m)", ylabel = "Y (m)",  guidefont=50)
    phi = range(0, stop = 2*pi, length = 100)
    x_obs = obs_info[1] .+ cos.(phi)
    y_obs = obs_info[2] .+ sin.(phi)
    plot!(x_obs, y_obs, line = (7, :black))
    x_obs2 = obs_info2[1] .+ cos.(phi)
    y_obs2 = obs_info2[2] .+ sin.(phi)
    plot!(x_obs2, y_obs2, line = (7, :black))
    savefig("test.png")
    
end
p = plot(size = [4000, 2000])
main()

