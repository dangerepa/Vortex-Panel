%% main_uvlm.m
% Aidan Bryce Sabin, May 2026
% Unsteady Vortex Panel Method, Aidan Bryce Sabin, May 2026
% Unsteady Vortex Lattice Method (2D)
% Test case: Sinusoidal pitching and/or plunging flat-plate or cambered airfoil
% -------------------------------------------------------------------------
clear; clc; close all;

%% Parameters 
% Flow
U_inf = 1.0; % Free-stream speed [m/s]
rho = 1.225; % Air density [kg/m^3]

% Geometry
c = 1.0; % Chord length [m]
N = 20; % Number of panels (bound vortices)

% Pitching motion: alpha(t) = alpha0 + alpha_amp * sin(omega*t + phi_p)
alpha0 = deg2rad(0); % Mean angle of attack [rad]
alpha_amp = deg2rad(5); % Pitch amplitude [rad]
% Want to test omega = 0.1, 0.4, 1
omega = 0.4; % Frequency [rad/s]
phi_p = pi; % Phase of pitching [rad]

% Plunging motion: z0(t) = h_amp * sin(omega*t)
h_amp = 0.1*c; % Plunge amplitude [m] (set 0 to disable)

% Pitch pivot location along chord (0 = LE, c = TE)
x_pivot = 0.25*c;

% Time
dt = 0.1; % Time step [s]
I_t_max = 400; % Number of time steps

% Derived
k = omega * c / (2*U_inf); % Reduced frequency
fprintf('Reduced frequency k = %.4f\n', k);

%% Airfoil definition 
% Approximates airfoil as camber line
USE_FLAT_PLATE = false;

if USE_FLAT_PLATE
    camber_x = [];
    camber_z = [];
else
    [x_af, y_af] = generate_naca4('0018', 100);
    [camber_x, camber_z] = load_airfoil([x_af, y_af], c);
end
%% Pre-allocate structures 
% We use struct arrays indexed by time snapshot.
% Fields are filled inside each function.
vb = struct('X', [], 'Z', [], 'Gamma', []);
vs = struct('X', [], 'Z', [], 'Gamma', []);
col = struct('X', [], 'Z', [], 'nhat', [], 'that', []);
cor = struct('X', [], 'Z', []);

% Storage for post-processing
t_hist = zeros(1, I_t_max);
L_hist = zeros(1, I_t_max);
M_hist = zeros(1, I_t_max);
CL_hist= zeros(1, I_t_max);

Gamma_wake_prev = 0; % Sum of all wake gammas at previous step (Kelvin)

%% Main time loop 
for I_t = 1:I_t_max
    t = I_t * dt;
    t_hist(I_t) = t;

 %  Step 1: Geometry (wing-fixed frame) 
    [cor, col_wf, vb_wf] = Geometry(c, N, I_t, cor, camber_x, camber_z);
 % col_wf / vb_wf hold collocation/bound-vortex data in wing-fixed frame

 %  Step 2: Kinematics 
    [x0, z0, alpha, x0dot, z0dot, alphadot] = ...
        Kinematics(t, dt, U_inf, alpha0, alpha_amp, omega, phi_p, h_amp);

 %  Step 3: Transform to inertial frame 
    [col, vb] = transform_coordinates(col_wf, vb_wf, x0, z0, alpha, ...
                                       x_pivot, I_t, col, vb);

 %  Step 4: Wake roll-up from PREVIOUS step (convect existing wake) 
 % Must happen before place_vortex so that vs(I_t) already contains the
 % convected old-wake positions when the new vortex is appended.
    if I_t > 1
        vs = wake_rollup(vb, vs, I_t-1, dt, U_inf);
 % vs(I_t).X/Z now holds convected positions for old wake vortices
    end

 %  Step 5: Place newly shed vortex 
 % TE_now and TE_prev computed purely from kinematics inside place_vortex.
    vs = place_vortex(vb, vs, U_inf, alpha0, alpha_amp, omega, phi_p, ...
                      h_amp, c, x_pivot, dt, I_t);

 %  Step 6: Build influence matrix A 
    A = Build_A(col_wf, vb_wf, vs, x0, z0, alpha, x_pivot, I_t, N);

 %  Step 7: Build RHS 
 % Gamma_wake_prev = sum of OLD wake vortices only (not the new one)
    RHS = Build_rhs(col_wf, vs, x0, z0, alpha, x0dot, z0dot, alphadot, ...
                    x_pivot, I_t, N, Gamma_wake_prev);

 %  Step 8: Solve linear system 
    Gamma_vec = Solve_system(A, RHS);

 % Store bound vortex strengths
    vb(I_t).Gamma = Gamma_vec(1:N);
 % Store newly shed vortex strength
    vs(I_t).Gamma(end) = Gamma_vec(N+1);

 % Update Kelvin sum: only OLD wake vortices for next step's RHS.
 % At next step, the vortex just shed becomes part of the old wake.
 % So Gamma_wake_prev for step I_t+1 = sum of all wake vortices at I_t.
    Gamma_wake_prev = sum(vs(I_t).Gamma);

 %  Step 9: Airloads 
    [L, M_LE, CL] = Calculate_airloads(col_wf, vb, vs, alpha, x0, z0, ...
                                        alphadot, x0dot, z0dot, x_pivot, ...
                                        I_t, dt, N, rho, U_inf, c);
    L_hist(I_t) = L;
    M_hist(I_t) = M_LE;
    CL_hist(I_t) = CL;

 % Progress
    if mod(I_t,20)==0
        fprintf('Step %4d / %4d  |  t = %.3f s  |  CL = %.4f\n', ...
                I_t, I_t_max, t, CL);
    end
end

%% Plot
plot_results(t_hist, L_hist, CL_hist, M_hist, vb, vs, I_t_max, ...
             U_inf, rho, c, omega);
