function BicycleModelDynamics(x, u)
    L_f = 1.0
    L_r = 1.0
    beta = atan(L_r / (L_r + L_f) * tan(u[2]))
    
    dx = [x[4] * cos(x[3] + beta),
            x[4] * sin(x[3] + beta),
            x[4] / L_r * sin(beta),
            u[1],
    ]
    return x + dx*dt
end

function KinematicBicycle(x, u)
    la = 1.0
    lb = 1.0
  
    dx= [(x[4]*cos(x[3] + (atan(la/(la+lb)*tan(u[2]))))),   # X position
    (x[4]*sin(x[3] + (atan(la/(la+lb)*tan(u[2]))))),   # Y position
    (x[4]*cos(atan(la/(la+lb)*tan(u[2])))/(la+lb)*tan(u[2])),      # Yaw Angle
    (u[1])]                                                  # Total Speed

    return dx
end

function DiscreteKinematicBycicle(x_t, u_t)
    A = ForwardDiff.jacobian((x) -> KinematicBicycle(x, u_t), x_t)
    B = ForwardDiff.jacobian((u) -> KinematicBicycle(x_t, u), u_t)
    Ad = exp(A*dt)
    f(v) = exp(A*v)
    Bd = quadgk(f, 0, dt)[1] * B
    return Ad*x_t + Bd*u_t
end