function [L, M_LE, CL] = Calculate_airloads(col_wf, vb, vs, alpha, ...
                           x0, z0, alphadot, x0dot, z0dot, x_pivot, ...
                           I_t, dt, N, rho, U_inf, c)
%CALCULATE_AIRLOADS  Compute lift and moment via unsteady Bernoulli.
% Aidan Bryce Sabin, May 2026
%  Katz & Plotkin (2001) Eq. 13.72 — pressure jump on panel j:
%
%    ΔP_j = ρ * [ V_tan_j * (Γ_j / Δl_j)  +  d/dt(Σ_{k=1}^{j} Γ_k) ]
%
%  where:
%    V_tan_j  = tangential velocity at panel j from KINEMATICS + WAKE only
%               (NOT from the bound vortex sheet — that is already in Γ_j)
%    Γ_j      = strength of bound vortex j
%    Δl_j     = panel arc length = c/N (flat plate)
%    d/dt(.)  = backward-difference of cumulative bound circulation
%
%  Lift per unit span (inertial z-component):
%    L = Σ_j  ΔP_j * Δl_j * n_{z,j,inertial}
%
%  Pitching moment about LE (positive nose-up):
%    M_LE = -Σ_j  ΔP_j * Δl_j * x_j * n_{z,j,inertial}
%
%  Inputs / outputs: see main_uvlm.m

ca = cos(alpha);  sa = sin(alpha);
R_ItoW = [ca, sa; -sa, ca]; % inertial -> wing-fixed
R_WtoI = [ca, -sa; sa, ca]; % wing-fixed -> inertial

dl_vec = col_wf.dl; % 1×N per-panel arc lengths

Gamma_bound = vb(I_t).Gamma; % 1×N, solved this step

%  Time derivative of cumulative bound circulation 
cum_now = cumsum(Gamma_bound);
if I_t > 1
    cum_prev = cumsum(vb(I_t-1).Gamma);
else
    cum_prev = zeros(1, N);
end
dGamma_dt = (cum_now - cum_prev) / dt; % 1×N  (added-mass term)

%  Kinematic velocity base in wing-fixed frame 
% x0dot = -U_inf  ->  -x0dot = +U_inf in wing frame (the free-stream)
% z0dot = plunge velocity
V_kin_base = R_ItoW * [-x0dot; -z0dot]; % includes U_inf via x0dot

%  All wake vortices (including newly shed) contribute to V_tan 
N_wake = length(vs(I_t).Gamma);

L = 0;
M_LE = 0;

for j = 1:N
    xj = col_wf.X(j);
    zj = col_wf.Z(j);
    tj = col_wf.that(:,j);
    nj_wf = col_wf.nhat(:,j);
    nj_I = R_WtoI * nj_wf; % normal in inertial frame

 %  Kinematic velocity at panel j (wing-fixed) 
 % Pitch contribution about x_pivot
    V_pitch = [-alphadot * zj;  alphadot * (xj - x_pivot)];
    V_kin = V_kin_base + V_pitch;

 %  Wake-induced velocity at panel j (wing-fixed) 
    u_w = 0;  w_w = 0;
    for jw = 1:N_wake
        p_I = [vs(I_t).X(jw) - x0;  vs(I_t).Z(jw) - z0];
        p_wf = R_ItoW * p_I;
        x_wv_wf = p_wf(1) + x_pivot;
        z_wv_wf = p_wf(2);
        [du, dw] = vor2d(vs(I_t).Gamma(jw), xj, zj, x_wv_wf, z_wv_wf);
        u_w = u_w + du;
        w_w = w_w + dw;
    end

 %  Tangential velocity (kinematics + wake; NOT bound vortex sheet) 
    V_tan = dot(V_kin + [u_w; w_w], tj);

 %  Pressure jump 
    dl = dl_vec(j);
    dP_j = rho * (V_tan * Gamma_bound(j) / dl + dGamma_dt(j));

 %  Lift and moment accumulation 
    L = L + dP_j * dl * nj_I(2);
    M_LE = M_LE - dP_j * dl * xj * nj_I(2);
end

%  Lift coefficient 
q_inf = 0.5 * rho * U_inf^2;
CL = L / (q_inf * c);
end
