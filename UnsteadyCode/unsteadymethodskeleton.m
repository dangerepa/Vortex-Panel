%% Background:
% Aidan Bryce Sabin, May 2026
% Want to save data as structures to make it easier to track time history, will do in a few steps. 
% From https://www.youtube.com/watch?v=4fCqD_lRd4Y
% Structures: 
%   vb: Bound Vortex Sheet Structure (inertial frame)
%   vs: shed vortices (inertial frame)
%   col: collocation points structure (inertial frame)
%   cor: corner point structure (wing-fixed frame)
% Example: vb(i).X(j) is i-th time snapshot of the x position of the j-th vortex 
% Skeleton of code: use several functions to be within a main script. 
%% Define α, U_\infty, rho, c (chord length), N (# of corner points), dt (time step), I_t_max (# of iterations)
%% 
for I_t = 1:I_t_max
    t = I_t*dt; 
    %Run steps 1 thru 8
end
%% Step 1
% Geometry.m
% Find Corner points, assign collocation points and bound vortices, n hat,
% t hat, in wing frame
% Same as steady method, but may need to store time history of these points
% Save to corner points structure, bound vortex structure, and collocation
% points structure. 
%% Step 2
% Kinematics.m
% Define wing kinematics, such as position, velocity, and angle of attack
% (angle between intertial frame and wing-fixed frame) and pitch rate. 
% Find: x0, z0, alpha, x0dot, z0dot, alpha dot
% Example: Impulsively started wing from rest (step input to freestream velocity)
%   x0 = -U_infty*t, alpha=cst. x0dot = -U_infty, z0=z0dot=alphadot=0
%% Step 3
% transform_coordinates.m
% transform all coordinates and vectors from wing-fixed frame to inertial
% save to col and vb structures. 
%% Step 4
% place_vortex.m
% This function defines the position of a newly shed vortex, but with an
% unknown strength (solved for later)
% Where to place vortex? Along TE path, 0.3*(Distance traveled by TE between t and t+Δt)
% save vortex position to vs (inertial position of vortex)
% every timestep, a new vortex element is shed -> vector within vs
% structure increases by 1 element every timestep. 
%% Step 5
% Build_A.m
% Similar to steady case
% Now we have to use :
% 1) Kelvin's circulation thm. dΓ/dt = 0 for all time. 
% 2) Find Γ of newly shed vortex
% From Kelvin, Γ1+Γ2+...Γ_N + Γ_new + Γ_wake(t-Δt)=0
% Essentially, all the gammas sum to the same gamma as the previous
% timestep. This will be represented by a row of 1s at the bottom of the A
% matrix. 
% The Matrix: first index is collocation pt, second is bound vortex
% [A][Γ]=[RHS]
%% Step 6
% Build_rhs.m
% Sets no flow-through BC. (applied at collocation pts in wing-fixed frame)
% Velocity contributions on each collocation point from 
%   a) Kinematics of wing
%   [u_k ; w_k] = R_{I->W}[-x0dot ; -z0dot] + [-αdot*col(i).z(j);αdot*col(i).x(j)]
%   " = Surge&Plunge + Pitching
%   b) Wake vortices
%   All wake vortices except newly shed vortex
%   Use vor2D.m function in the wing-fixed frame
%   Example: velocity contribution of wake vortex j on collocation point i
%   in wing-fixed frame
% [a_{w{j}, w_{w_j}] = vor2d(Γ_{w,j},x_i,z_i,w_{w,j},z_{w,j})
% The velocity contribution at collocation point i by all wake vortices
% [u_w,w_w] = [sum_{j=1}^{N_w} u_{w_j}, sum_{j=1}^{N_w} w_{w_j}]
% Do at each collocation point
% For collocation point i, RHS_i = -[u_k+u_w,w_k+w_w]\cdot n_i
% where n_i is in the wing-fixed frame
% and RHS_{N_c+1} = -Γ_w(t_Δt)
%% Step 7
% Solve_system.m
% Gamma = [Γ_1;Γ_2;...;Γ_{Nv};Γ_{new}] = A\RHS
%% Step 8
% Calculate_airloads.m
% Require ΔP across camber line: Use unsteady Bernoulli. 
% ΔP_j = ρ(([u_k+u_w, w_k+w_w]\cdot t_j)Γ_j/Δl_j+\partial/{\partial t}(\sum_{k=1}^j Γ_k))
% Interpretation: ΔP for panel j or vortex j, t is tangent unit vector at
% each panel j or bound vortex j. l_j is length of panel b/wn 2 adjacent
% vortices. Γ_j is strength of vortex j, and Γ_k is vortex strength of kth
% vortex. As j increases, you sum more and more vortices. 
%% Lift
% Lift is found by numericall integrating integrating pressure difference
% and obtaining the vertical component in inertial frame
% L_j = ΔP_j*Δl_j*n_{z_j,inertial}, L = sum(L_j's)
% L_A = sum_{j=1}^{Nv} ρ([u_k+u_w, w_k+w_w]\cdot t_j)Γ)j*n_{z_{j,inertial}}
% L_B = sum_{j=1}^{Nv} ρ*Δl_j*\partial/{\partial t} * sum_{k=1}^j Γ_k*n_{z_{j,inertial}}
% L = L_A+L_B
%% Pitching moment about LE
%M_{LE} = -\sum_{j=1}^{Nv}\Delta P_{j}\Delta l_{j}x_{j}n_{z_{j,inertial}}
% where x is in the wing-fixed frame
%% Step 9
% wake_rollup.m
% determine motion of wake (velocity), update wake position, done in
% inertial frame. Done in 2 parts
% A) For kth wake vortex, determine velocity induced at its position by 
%   i) Nv bound vortices
%   ii) other vortices in the wake. (include gamma of newly shed vortex)
% B) Convecting vortices in wake (Euler's method)
% dx/dt = \frac{x^{n+1}-x^{n}}{\Delta t}
% vs(i+1).X(k) = vs(i).X(k)+U_k*Δt
% vs(i+1).Y(k) = vs(i).Y(k)+w_k*Δt


% Vor2d: 
% SUBROUTINE VOR2D(X,Z,XI,Zl,GAMMA,U,W)
% CALCULATES INFLUENCE OF VORTEX AT (XI,Zl)
% PAY=3.141592654
% U=0.0
% W=0.0
% RX=X-X1
% RZ=Z-Z1
% R=SQRT(RX**2+RZ**2)
% IF(R.LT.0.001)GOTO 1
% V=0.5/PAY*GAMMA/R
% U=V*(RZ/R)
% W=V*(-RX/R)
% CONTINUE
% RETURN
% END