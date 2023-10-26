function BicycleModelDynamics(x, u)
   #=
    L_f = 1.0
    L_r = 1.0
    beta = atan(L_r / (L_r + L_f) * tan(u[2]))
    
    dx = [x[4] * cos(x[3] + beta),
            x[4] * sin(x[3] + beta),
            x[4] / L_r * sin(beta),
            u[1],
    ]
    return x + dx*dt
   =#
    
    la = 1.0;
    lb = 1.0;
    x_pos =(x[4]*cos(x[3] + (atan(la/(la+lb)*tan(u[2])))))   # X position
    y_pos =(x[4]*sin(x[3] + (atan(la/(la+lb)*tan(u[2])))))   # Y position
    yaw =(x[4]*cos(atan(la/(la+lb)*tan(u[2])))/(la+lb)*tan(u[2])   )      # Yaw Angle
    accel =u[1]                                                  # Total Speed
    dx = [x_pos, y_pos, yaw, accel]
    return x+dx*dt
    
end