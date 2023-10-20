using LinearAlgebra
using JuMP
using Ipopt
using Mosek
using MosekTools
using COPT
using SCS
using ForwardDiff

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

function CovarianceControl(x_ref, u_ref)
    A = zeros(Float64, (N, x_dim, x_dim))
    B = zeros(Float64, (N, x_dim, u_dim))

    for k in 1:N
        local Ak = ForwardDiff.jacobian((x) -> BicycleModelDynamics(x, u_ref[k,:]), x_ref[k,:])
        local Bk = ForwardDiff.jacobian((u) -> BicycleModelDynamics(x_ref[k,:], u), u_ref[k,:])
        A[k,:,:] = Ak
        B[k,:,:] = Bk
    end

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
    print(Beta)
    # model = Model(Ipopt.Optimizer)
    # model = Model(Mosek.Optimizer)
    # model = Model(COPT.ConeOptimizer)
    model = Model(SCS.Optimizer)
    # set_silent(model)
    @variable(model, K[1:N*u_dim, 1:(N+1)*x_dim])
    # @constraint(model, En *(I+B*K)*B*sigma_control_bar*B'*(I+B*K)'*En' .<= sigma_xf)
    @constraint(model, [sigma_xf En*(I+Beta*K)*Beta*sqrt(sigma_control_bar); sqrt(sigma_control_bar)*Beta'*(I+Beta*K)'*En' I] in PSDCone())
    obj = tr(((I + Beta*K)' * Q_bar * (I + Beta*K) + K'R_bar*K)*Beta*sigma_control_bar*Beta')
    @objective(model, Min, obj)
    optimize!(model)

    return A, B, value.(K)
end