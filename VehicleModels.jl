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