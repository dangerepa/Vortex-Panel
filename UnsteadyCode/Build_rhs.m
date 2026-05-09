function RHS = Build_rhs(col_wf, vs, x0, z0, alpha, x0dot, z0dot, ...
                          alphadot, x_pivot, I_t, N, Gamma_wake_prev)
%BUILD_RHS  Assemble the right-hand side vector.
% Aidan Bryce Sabin, May 2026
%  RHS_i  = -( V_kin + V_wake ) · n_i      for i = 1..N  (BC rows)
%  RHS_N+1 = -Γ_wake(t-Δt)                               (Kelvin row)
%
%  All velocities are evaluated in the wing-fixed frame.
%
%  Kinematic velocity at collocation point i (wing-fixed frame):
%    V_kin = R_{I->W} * [-x0dot; -z0dot]          (surge + plunge)
%          + [-alphadot*z_i_wf; alphadot*x_i_wf]  (pitch)
%
%  Wake velocity: sum over ALL wake vortices except the newly shed one.
%
%  Inputs:
%    col_wf           - collocation points in wing-fixed frame
%    vs               - wake vortex structure (inertial)
%    x0, z0           - pivot position (inertial)
%    alpha            - pitch angle [rad]
%    x0dot, z0dot     - pivot velocity (inertial)
%    alphadot         - pitch rate [rad/s]
%    x_pivot          - pivot location along chord (wing-fixed)
%    I_t              - time index
%    N                - number of panels
%    Gamma_wake_prev  - total wake circulation at previous time step (Kelvin)
%
%  Outputs:
%    RHS  - (N+1) × 1 right-hand side vector

RHS = zeros(N+1, 1);

ca = cos(alpha);  sa = sin(alpha);
R_ItoW = [ca, sa; -sa, ca]; % Rotation: inertial -> wing-fixed

% Velocity of pivot in inertial frame, brought to wing-fixed frame
% The wing sees a free-stream contribution = -V_pivot in wing frame
V_pivot_I  = [x0dot; z0dot];
V_kinematic_base = R_ItoW * (-V_pivot_I); % surge + plunge in wing frame

% Number of OLD wake vortices (exclude newly shed vortex at index end)
N_wake = length(vs(I_t).Gamma) - 1; % newly shed is the last entry

for i = 1:N
    xi_wf = col_wf.X(i);
    zi_wf = col_wf.Z(i);
    ni    = col_wf.nhat(:,i);

 %  Kinematic velocity at this collocation point 
 % Pitch contribution: rotation about pivot in wing frame
    V_pitch_wf = [-alphadot * zi_wf; alphadot * (xi_wf - x_pivot)];

    V_kin = V_kinematic_base + V_pitch_wf;

 %  Wake velocity at this collocation point 
    u_wake = 0;  w_wake = 0;
    for jw = 1:N_wake
 % Wake vortex position in inertial frame
        x_wv_I = vs(I_t).X(jw);
        z_wv_I = vs(I_t).Z(jw);
        Gam_wv = vs(I_t).Gamma(jw);

 % Transform wake vortex to wing-fixed frame
        p_I  = [x_wv_I - x0; z_wv_I - z0];
        p_wf = R_ItoW * p_I;
        x_wv_wf = p_wf(1) + x_pivot;
        z_wv_wf = p_wf(2);

        [du, dw] = vor2d(Gam_wv, xi_wf, zi_wf, x_wv_wf, z_wv_wf);
        u_wake = u_wake + du;
        w_wake = w_wake + dw;
    end

    V_wake = [u_wake; w_wake];

 %  RHS for this collocation point 
    RHS(i) = -dot(V_kin + V_wake, ni);
end

%  Kelvin row
RHS(N+1) = -Gamma_wake_prev;
end
