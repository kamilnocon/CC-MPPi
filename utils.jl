using LinearAlgebra
using JuMP
using Ipopt
using Mosek
using MosekTools
#using COPT
using SCS
using ForwardDiff
using QuadGK
using DelimitedFiles

function cost(x0, u0, obs_info, obs_info2, eps, R)
    state_cost = 0
    if ((x0[1] < -30.0) | (x0[1] > 30.0)) & ((x0[2] < -5) | (x0[2] > 5))
        state_cost = state_cost + 200000.0
    end
    
    # TODO: Sigmoid or Tanh function for obstacles + boundaries
    if sqrt((x0[1]-obs_info[1])^2+(x0[2]-obs_info[2])^2) < obs_info[3]+0.1
        state_cost = state_cost + 1000.0*712.5
    end
    if sqrt((x0[1]-obs_info2[1])^2+(x0[2]-obs_info2[2])^2) < obs_info2[3]+0.1
        state_cost = state_cost + 1000.0*712.5
    end
    
    #state_cost += 150000*max(0, obs_info[3]-sqrt((x0[1]-obs_info[1])^2+(x0[2]-obs_info[2])^2))
    #state_cost += 150000*max(0, obs_info2[3]-sqrt((x0[1]-obs_info2[1])^2+(x0[2]-obs_info2[2])^2))
    
    if (x0[3] < -2*pi) | (x0[3] > 2*pi)
        state_cost += 7125
    end
    if (u0[1] < -4) | (u0[1] > 4)
        state_cost += 7125
    end
    if (x0[5] < -0.62) | (x0[5] > 0.62)
        state_cost += 7125
    end
    if (x0[4] < 2) | (x0[4] > 10)
        state_cost += 7125
    end
    if (u0[2] < -0.56) | (u0[2] > 0.56)
        state_cost += 7125
    end
    
    state_cost += 50.0*(u0[1])^2
    state_cost += 50.0*(u0[2])^2
    state_cost += 1800.0*(x0[2]-0.0)^2
    state_cost += 500.0*(x0[1]-11.0)^2
    
    state_cost += 1000/x0[4]

    state_cost += 100*abs(0.0-x0[2])

    sampling_cost = (((1-nu^-1)/2)*transpose(eps)*R*eps)[1]+(transpose(u0)*R*eps)[1] + (1/2)*(transpose(u0)*R*u0)[1]
    rollout_cost = state_cost + sampling_cost
    return rollout_cost
end

function state_cost(x0, u0)
    state_cost = 0
    if ((x0[1] < -30.0) | (x0[1] > 30.0)) & ((x0[2] < -5) | (x0[2] > 5))
        state_cost = state_cost + 200000.0
    end
    
    if sqrt((x0[1]-obs_info[1])^2+(x0[2]-obs_info[2])^2) < obs_info[3]+0.1
        state_cost = state_cost + 1000.0*712.5
    end
    if sqrt((x0[1]-obs_info2[1])^2+(x0[2]-obs_info2[2])^2) < obs_info2[3]+0.1
        state_cost = state_cost + 1000.0*712.5
    end
    
    #state_cost += 150000*max(0, obs_info[3]-sqrt((x0[1]-obs_info[1])^2+(x0[2]-obs_info[2])^2))
    #state_cost += 150000*max(0, obs_info2[3]-sqrt((x0[1]-obs_info2[1])^2+(x0[2]-obs_info2[2])^2))
    
    if (x0[3] < -2*pi) | (x0[3] > 2*pi)
        state_cost += 7125
    end
    if (u0[1] < -4) | (u0[1] > 4)
        state_cost += 7125
    end
    if (x0[5] < -0.62) | (x0[5] > 0.62)
        state_cost += 7125
    end
    if (x0[4] < 2) | (x0[4] > 10)
        state_cost += 7125
    end
    if (u0[2] < -0.56) | (u0[2] > 0.56)
        state_cost += 7125
    end
    
    state_cost += 50.0*(u0[1])^2
    state_cost += 50.0*(u0[2])^2
    state_cost += 1800.0*(x0[2]-0.0)^2
    state_cost += 500.0*(x0[1]-11.0)^2
    
    state_cost += 1000/x0[4]
    return state_cost
end

function terminal_cost(x0, x)
    s = abs(x0[1]-x)
    e = abs(0.0-x0[2])
    term_cost = 3.3*(1-s)+4000*e^2
    return term_cost
end

function optimal_control(Sm_list, Um_list, lambda, M)
    b = minimum(Sm_list)
    V = zero(Um_list[1,:,:])
    sum = 0
    for m in 1:M
        wm = exp((-1/lambda)*(Sm_list[m]-b))
        V += wm.*Um_list[m, :, :]
        sum += wm
    end
    V = V ./ sum
    return V
end

function CovarianceControl(x_ref, u_ref, opt_model)
    A = zeros(Float64, (N, x_dim, x_dim))
    B = zeros(Float64, (N, x_dim, u_dim))
    for k in 1:N
        # local Ak = ForwardDiff.jacobian((x) -> KinematicBicycle(x, u_ref[k,:]), x_ref[k,:])
        # local Bk = ForwardDiff.jacobian((u) -> KinematicBicycle(x_ref[k,:], u), u_ref[k,:])
        # A[k,:,:] = exp(Ak*dt)
        # f(v) = exp(Ak*v)
        # B[k,:,:] = quadgk(f, 0, dt)[1] * Bk
        
        local Ak = ForwardDiff.jacobian((x) -> BicycleModelDynamics(x, u_ref[k,:]), x_ref[k,:])
        local Bk = ForwardDiff.jacobian((u) -> BicycleModelDynamics(x_ref[k,:], u), u_ref[k,:])
        A[k,:,:] = Ak
        B[k,:,:] = Bk

        local Qk = ForwardDiff.hessian((x) -> state_cost(x, u_ref[k,:]), x_ref[k,:])
        local Rk = ForwardDiff.hessian((u) -> state_cost(x_ref[k,:], u), u_ref[k,:])
        for k in 1:N
            Q_bar[(k-1)*x_dim + 1:k*x_dim , (k-1)*x_dim + 1:k*x_dim] = Qk * 0.0001
            R_bar[(k-1)*u_dim + 1:k*u_dim , (k-1)*u_dim + 1:k*u_dim] = Rk * 0.0001
        end
    end
    Q_bar[N*x_dim + 1:(N+1)*x_dim , N*x_dim + 1:(N+1)*x_dim] = ForwardDiff.hessian((x) -> terminal_cost(x, x_ref[1,1]), x_ref[N+1,:]) 

    Beta = zeros(Float64, ((N+1)*x_dim, (N)*u_dim))
    for k in 1:N
        for j in 1:k
            local Bk1k0 = B[j,:,:]
            for z in j:k
                Bk1k0 = A[z,:,:] * Bk1k0
            end
            Beta[k*x_dim+1:(k+1)*x_dim, (j-1)*u_dim+1:(j)*u_dim] = Bk1k0
        end
    end
    #set_silent(model)
    
    #@constraint(model, En*(I+Beta*K)*Beta*sigma_control_bar*Beta'*(I+Beta*K)'*En' .<= sigma_xf)
    #@constraint(model, [sigma_xf En*(I+Beta*K)*Beta*sqrt(sigma_control_bar); sqrt(sigma_control_bar)*Beta'*(I+Beta*K)'*En' I] in PSDCone())
    @constraint(opt_model, [sigma_xf En*(I+Beta*K)*Beta*sqrt(sigma_control_bar); sqrt(sigma_control_bar)*Beta'*(I+Beta*K)'*En' I] >= 0, PSDCone())
    obj = tr(((I + Beta*K)' * Q_bar * (I + Beta*K) + K'R_bar*K)*Beta*sigma_control_bar*Beta')
    @objective(opt_model, Min, obj)
    optimize!(opt_model)
    for i in 1:N
        for j in 1:N+1
            if (i+1 < j) | (i+1 > j)
                value.(K)[(i-1)*u_dim+1:i*u_dim, (j-1)*x_dim+1:j*x_dim] .= 0
            end
        end
    end
    writedlm("K.csv", value.(K), ',')
    return A, B, value.(K)
end