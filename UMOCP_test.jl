using UMOCP
using Test

@testset "UMOCP.jl" begin

    XL = [-2.6, 0, NaN, NaN, -6.283185307179586, -0.62, 10, -2.6]
    XU = [1.801, 300, NaN, NaN, 6.283185307179586, 0.62, 20,  2.6]
    CL=[-2.1, -6.2]
    CU=[1.9, 6.2]
    n = define(numStates=8,numControls=2,X0=[1.8, 0., 0.0, 0.0, pi/2, 0.0, 15., 0.],XF=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], XL=XL, XU=XU, CL=CL,CU=CU);

    goal = [1.8, 200, 5]
    # goal = [-1.8, 100, 5] # if you want to move to left lane by the end
    obs_info = [1.8 50.7368 0.3 0.9 0.0 5.5];
    obs_num = size(obs_info, 1)

    states!(n,[:x,:y,:v,:r,:psi,:sa,:ux,:ax];descriptions=["x(t)","y(t)","v(t)","r(t)","psi(t)","sa(t)","ux(t)","ax(t)"]);
    controls!(n,[:jx,:sr];descriptions=["jx(t)","sr(t)"]);
    dx=ThreeDOFBicycle_expr(n);
    dynamics!(n,dx);
    configure!(n,N=60;(:integrationScheme=>:bkwEuler),(:finalTimeDV=>true),(:solverSettings => (:name => :Ipopt)));
    # configure!(n,N=60;(:integrationScheme=>:bkwEuler),(:finalTimeDV=>false), (:tf => 6));
    # configure!(n;(:Nck=>[60]),(:integrationScheme=>:lgrImplicit),(:finalTimeDV=>true),(:solverSettings => (:name => :Ipopt)));
    # configure!(n,N=60;(:integrationScheme=>:mpcol),(:finalTimeDV=>false), (:tf => 6));
    # configure!(n,N=40;(:integrationScheme=>:mpcol),(:finalTimeDV=>true), (:solverSettings => (:name => :MadNLP)));
    #

    n.s.ocp.evalConstraints = true
    lsm = 4.66 #to leave 2m gap between the tip of the vehicle and the behind of the bicycle
    ssm = 3.3 #to move the vehicle to the center of the left lane.
    x = n.r.ocp.x[:,1];y = n.r.ocp.x[:,2];ux = n.r.ocp.x[:,7];psi = n.r.ocp.x[:,5];# pointers to JuMP variables
    timeSeq = n.ocp.tV[:,1];

    #obs_con = @NLconstraint(n.ocp.mdl, [i=1:n.ocp.state.pts-1], 1 <= ((x[(i+1)]-timeSeq[i+1]*c["obstacle"]["vx"][1]-c["obstacle"]["x0"][1])^2)/((c["obstacle"]["saxis"][1]+ssm)^2) + ((y[(i+1)]-timeSeq[i+1]*c["obstacle"]["vy"][1]-c["obstacle"]["y0"][1])^2)/((c["obstacle"]["laxis"][1]+lsm)^2))

    obs_ind = 1
    obs_con1 = @NLconstraint(n.ocp.mdl, [i=1:n.ocp.state.pts-1], 1 <= ((x[(i+1)]-timeSeq[i+1]*obs_info[obs_ind, 5]-obs_info[obs_ind, 1])^2)/((obs_info[obs_ind, 3]+ssm)^2) + ((y[(i+1)]-timeSeq[i+1]*obs_info[obs_ind, 6]-obs_info[obs_ind, 2])^2)/((obs_info[obs_ind, 4]+lsm)^2))
    newConstraint!(n,obs_con1,:obs_con1)


    ux_con2 = @NLconstraint(n.ocp.mdl, [i=1:n.ocp.state.pts-1], (ux[i+1]-15)^2<=0.01)
    newConstraint!(n,ux_con2,:ux_con2)

    # goal_con = @NLconstraint(n.ocp.mdl, [i=1:n.ocp.state.pts-1], (x[n.ocp.state.pts]-goal[1])^2+(y[n.ocp.state.pts]-goal[2])^2<=goal[3]^2)
    # goal_con = @NLconstraint(n.ocp.mdl, [i=1:n.ocp.state.pts-1], (y[n.ocp.state.pts]-goal[2])>=200)
    goal_con = @NLconstraint(n.ocp.mdl, [i=n.ocp.state.pts-5:n.ocp.state.pts-1], ((psi[i]-pi/2)^2<=0.001))
    newConstraint!(n,goal_con,:goal_con)


    obj = integrate!(n,:( 10*sr[j]^2+ 5*(x[j]-1.8)^2 + 10*v[j]^2));
    print("here\n")

    @NLobjective(n.ocp.mdl, Min, obj+n.ocp.tf+(n.r.ocp.x[end,1]-goal[1])^2 + (n.r.ocp.x[end,2]-goal[2])^2)
    NLoptimize!(n);
end