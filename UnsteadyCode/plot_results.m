function plot_results(t_hist, L_hist, CL_hist, M_hist, vb, vs, I_t_max, ...
                      U_inf, rho, c, omega)
%PLOT_RESULTS  Visualise UVLM outputs.
% Aidan Bryce Sabin, May 2026
%  Produces four figures:
%    1. CL vs non-dimensional time (tU/c)
%    2. Lift L vs time
%    3. Pitching moment about LE vs time
%    4. Final wake shape (vortex positions)

tUc = t_hist * U_inf / c; % non-dimensional time

%%  Figure 1: CL vs tU/c 
figure(1); clf;
plot(tUc, CL_hist, 'b-', 'LineWidth', 1.5);
xlabel('tU_\infty/c');
ylabel('C_L');
title('Lift Coefficient vs Non-dimensional Time');
grid on;

%%  Figure 2: Lift vs time 
figure(2); clf;
plot(t_hist, L_hist, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Lift per unit span (N/m)');
title('Lift vs Time');
grid on;

%%  Figure 3: Pitching moment vs time 
figure(3); clf;
plot(t_hist, M_hist, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('M_{LE} (N·m/m)');
title('Pitching Moment about LE vs Time');
grid on;

%%  Figure 4: Wake geometry at final time step 
figure(4); clf;
hold on; box on;

% Bound vortices
plot(vb(I_t_max).X, vb(I_t_max).Z, 'bs-', 'MarkerFaceColor', 'b', ...
     'DisplayName', 'Bound vortices');

% Wake vortices
if ~isempty(vs(I_t_max).X)
    scatter(vs(I_t_max).X, vs(I_t_max).Z, 20, vs(I_t_max).Gamma, 'filled', ...
            'DisplayName', 'Wake vortices');
    colorbar; colormap(jet);
    clim_val = max(abs(vs(I_t_max).Gamma));
    if clim_val > 0
        clim([-clim_val, clim_val]);
    end
end

xlabel('x (m)');
ylabel('z (m)');
title('Wake Shape at Final Time Step');
legend('Location','best');
axis equal;
grid on;

end
