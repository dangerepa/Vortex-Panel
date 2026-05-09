function vs = place_vortex(vb, vs, U_inf, alpha0, alpha_amp, omega, phi_p, ...
                            h_amp, c, x_pivot, dt, I_t)
%PLACE_VORTEX  Shed a new vortex from the trailing edge.
% Aidan Bryce Sabin, May 2026
%  The newly shed vortex is placed at 30% of one dt of TE travel downstream
%  of the current TE position (Katz & Plotkin recommendation):
%
%    P_new = TE_now + 0.3 * (TE_now - TE_prev)
%
%  Both TE_now and TE_prev are computed PURELY from body kinematics:
%
%    TE(t) = [x0(t); z0(t)] + R_{W->I}(alpha(t)) * [c - x_pivot; 0]
%
%  where R_{W->I} = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]
%
%  This is the correct approach. Using vs(I_t-1).X(end) as TE_prev is wrong
%  because that point is already offset from the TE by the 0.3-fraction rule.
%
%  Inputs:
%    vb                          - bound vortex structure (unused, kept for call signature)
%    vs                          - wake vortex structure
%    U_inf, alpha0, alpha_amp,
%    omega, phi_p, h_amp         - kinematic parameters
%    c                           - chord length
%    x_pivot                     - pivot x-location in wing-fixed frame
%    dt                          - time step
%    I_t                         - current time index
%
%  Outputs:
%    vs  - updated; old-wake positions already set by wake_rollup;
%          new vortex appended at end with Gamma = 0.

t_now  = I_t * dt;
t_prev = (I_t - 1) * dt;

%  TE in inertial frame at t_now and t_prev 
TE_now = TE_inertial(t_now,  U_inf, alpha0, alpha_amp, omega, phi_p, h_amp, c, x_pivot);
TE_prev = TE_inertial(t_prev, U_inf, alpha0, alpha_amp, omega, phi_p, h_amp, c, x_pivot);

%  New vortex: 30% of one TE-travel step ahead of current TE 
frac  = 0.3;
TE_step = TE_now - TE_prev;
x_new = TE_now(1) - frac * TE_step(1);
z_new = TE_now(2) - frac * TE_step(2);
%  Append to wake snapshot I_t 
% Convected positions 1..I_t-1 are already in vs(I_t) from wake_rollup.
if I_t == 1
    vs(I_t).X = x_new;
    vs(I_t).Z = z_new;
    vs(I_t).Gamma = 0;
else
    n_old = length(vs(I_t-1).Gamma);
    vs(I_t).Gamma(1:n_old) = vs(I_t-1).Gamma(1:n_old);
    vs(I_t).X(n_old+1) = x_new;
    vs(I_t).Z(n_old+1) = z_new;
    vs(I_t).Gamma(n_old+1) = 0;
end
end

% =========================================================================
function TE = TE_inertial(t, U_inf, alpha0, alpha_amp, omega, phi_p, h_amp, c, x_pivot)
%TE_INERTIAL  Trailing-edge position in inertial frame at time t.
%
%  In the wing-fixed frame the TE sits at (c - x_pivot, 0) relative to the
%  pivot.  The rotation R_{W->I} and translation (x0, z0) map it to inertial:
%
%    TE = [x0; z0] + R_{W->I} * [c - x_pivot; 0]
%
%  R_{W->I} = [ cos(alpha), -sin(alpha) ]
%             [ sin(alpha), cos(alpha) ]
%
%  So:  TE_x = x0 + (c - x_pivot)*cos(alpha)
%       TE_z = z0 + (c - x_pivot)*sin(alpha)
%
%  Sign note: positive alpha (nose up) rotates the TE downward in z only
%  if we use a different convention.  Here alpha is the angle of the wing
%  x-axis relative to the inertial x-axis, so the standard rotation applies
%  and TE_z = z0 + (c-x_pivot)*sin(alpha) (positive alpha -> TE moves up).

x0 = -U_inf * t;
z0 =  h_amp * sin(omega * t);
alpha = alpha0 + alpha_amp * sin(omega * t + phi_p);

d = c - x_pivot; % distance from pivot to TE along wing x-axis
TE = [x0 + d * cos(alpha);
      z0 + d * sin(alpha)];
end
