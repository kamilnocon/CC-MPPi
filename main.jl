using UMOCP
using JuMP
using Ipopt
using Plots
include("parameters.jl")
include("VehicleModels.jl")

obs_pos = [20, 0, 3.5]
x0 = 0; y0 = 1e-3; v0 = 0; r0 = 0; psi0 = 0;
init_states = [x0, y0, v0, r0, psi0, x0, y0, v0, r0, psi0, 15.0]
XL = [NaN, -0.5, -5, -pi/3, -pi/3, NaN, -0.5, -5, -pi/3, -pi/3, 5]
XU = [NaN, NaN, 5, pi/3, pi/3, NaN, NaN, 5, pi/3, pi/3, 25]
CL = [-pi/12, -pi/12]
CU = [pi/12, pi/12]

X0 = init_states
n = define(numStates=11,numControls=2,X0=X0,XF=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], XL=XL, XU=XU, CL=CL,CU=CU);

states!(n,[:x1,:y1,:v1,:r1,:psi1,:x2,:y2,:v2,:r2,:psi2,:ux];descriptions=["x1(t)","y1(t)","v1(t)","r1(t)","psi1(t)","x2(t)","y2(t)","v2(t)","r2(t)","psi2(t)","ux(t)"]);
controls!(n, [:sa1, :sa2]; descriptions=["sa1(t)", "sa2(t)"])
dx = ContingentVehicleModel(n)
dynamics!(n, dx)

NckPoint = 20;
PredictionTimeHorizon = 3;
configure!(n,N=NckPoint;(:integrationScheme=>:bkwEuler),(:tf=>PredictionTimeHorizon), (:solverSettings => (:name => :Ipopt)));


x1 = n.r.ocp.x[:,1];y1 = n.r.ocp.x[:,2];x2 = n.r.ocp.x[:, 6]; y2 = n.r.ocp.x[:, 7]; sa1 = n.r.ocp.u[:,1];sa2 = n.r.ocp.u[:,2]# pointers to JuMP variables


contingentDataPoints = 5;
CTRLCon = @NLconstraint(n.ocp.mdl, [j = 1:contingentDataPoints], sa1[j] == sa2[j])
newConstraint!(n,CTRLCon,:CTRLCon)
obs1 = @NLconstraint(n.ocp.mdl, [j = 1:n.ocp.state.pts], ((x1[j] - obs_pos[1])^2 / obs_pos[3]^2 + (y1[j] - obs_pos[2])^2 / obs_pos[3]^2) >= 1 )
newConstraint!(n,obs1,:obs1)
obs2 = @NLconstraint(n.ocp.mdl, [j = 1:n.ocp.state.pts], ((x2[j] - obs_pos[1])^2 / obs_pos[3]^2 + (y2[j] - obs_pos[2])^2 / obs_pos[3]^2) >= 1 )
newConstraint!(n,obs2,:obs2)


obj = integrate!(n,:( 10* (0.5 * sa1[j]^2 + 0.5 * sa2[j]^2 ) + 1 * log((y1[j]^2 + y2[j]^2 + 0.01))));#+1.0*y[j]^2
@NLobjective(n.ocp.mdl, Min, obj)

NLoptimize!(n)

p = plot(layout=(1,2), size = [4000, 2000])

plot!(subplot = 1, n.r.ocp.X[:, 1], n.r.ocp.X[:, 2], tickfontsize=30, line = (7, :blue),     xlabel = "X (m)", ylabel = "Y (m)",  guidefont=50       )
plot!(subplot = 1, n.r.ocp.X[:, 6], n.r.ocp.X[:, 7], tickfontsize=30, line = (7, :red),     xlabel = "X (m)", ylabel = "Y (m)",  guidefont=50       )

plot!(subplot = 2, n.r.ocp.tctr[:, 1], n.r.ocp.U[:, 1], tickfontsize=30, line = (7, :blue),     xlabel = "t (s)", ylabel = "SA (rad)",  guidefont=50       )
plot!(subplot = 2, n.r.ocp.tctr[:, 1], n.r.ocp.U[:, 2], tickfontsize=30, line = (7, :red),     xlabel = "t (s)", ylabel = "SA (rad))",  guidefont=50       )