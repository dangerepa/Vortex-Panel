function vs = wake_rollup(vb, vs, I_t, dt, U_inf)
% WAKE_ROLLUP  Convect all wake vortices to their next-step positions.
% Aidan Bryce Sabin, May 2026
%
%  Uses Euler's method (first-order explicit time integration):
%    x^{n+1} = x^n + U_k * Δt
%    z^{n+1} = z^n + W_k * Δt
%
%  The velocity at each wake vortex position is the superposition of:
%    (i)  All N bound vortices
%    (ii) All other wake vortices (including newly shed)
%
%  Results are stored in vs(I_t+1) to set up the next time step.
%
%  Inputs:
%    vb   - bound vortex structure, snapshot I_t
%    vs   - wake vortex structure, snapshot I_t (all Gammas solved)
%    I_t  - current time index
%    dt   - time step
%    U_inf - free-stream speed (not used directly; already encoded in kinematics)
%
%  Outputs:
%    vs   - updated; vs(I_t+1).X/Z set to convected positions

N_wake = length(vs(I_t).Gamma);
N_bv   = length(vb(I_t).Gamma);

%  Compute induced velocity at each wake vortex position 
U_k = zeros(1, N_wake);
W_k = zeros(1, N_wake);

for k = 1:N_wake
    xk = vs(I_t).X(k);
    zk = vs(I_t).Z(k);

    u_total = 0;
    w_total = 0;

  % Contribution from all bound vortices
    for j = 1:N_bv
        [du, dw] = vor2d(vb(I_t).Gamma(j), xk, zk, ...
                         vb(I_t).X(j), vb(I_t).Z(j));
        u_total = u_total + du;
        w_total = w_total + dw;
    end

  % Contribution from all other wake vortices
    for jw = 1:N_wake
        if jw == k, continue; end % skip self
        [du, dw] = vor2d(vs(I_t).Gamma(jw), xk, zk, ...
                         vs(I_t).X(jw), vs(I_t).Z(jw));
        u_total = u_total + du;
        w_total = w_total + dw;
    end

    U_k(k) = u_total;
    W_k(k) = w_total;
end

%  Euler step
vs(I_t+1).X     = vs(I_t).X     + U_k * dt;
vs(I_t+1).Z     = vs(I_t).Z     + W_k * dt;
vs(I_t+1).Gamma = vs(I_t).Gamma; % strengths don't change (frozen wake)
end
