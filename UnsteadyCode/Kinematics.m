function [x0, z0, alpha, x0dot, z0dot, alphadot] = ...
    Kinematics(t, dt, U_inf, alpha0, alpha_amp, omega, phi_p, h_amp)
%KINEMATICS  Define wing position, velocity, and attitude.
% Aidan Bryce Sabin, May 2026
%  Sinusoidal pitching + plunging about the pivot point.
%
%  Motion definitions (inertial frame, pivot as reference):
%    Plunge:  z0(t)     =  h_amp  * sin(omega*t)
%    Pitch:   alpha(t)  =  alpha0 + alpha_amp * sin(omega*t + phi_p)
%    Surge:   x0(t)     = -U_inf * t   (wing translates into free-stream)
%
%  Sign convention:
%    x0, z0   - position of the pivot/reference point in inertial frame
%    alpha     - angle from inertial x-axis to wing x_w-axis (positive nose-up)
%    dots      - time derivatives of the above
%
%  Inputs:
%    t        - current time [s]
%    dt       - time step (unused here, kept for finite-difference use)
%    U_inf    - free-stream speed
%    alpha0   - mean pitch angle [rad]
%    alpha_amp- pitch amplitude [rad]
%    omega    - angular frequency [rad/s]
%    phi_p    - pitch phase offset [rad]
%    h_amp    - plunge amplitude [m]
%
%  Outputs:
%    x0, z0        - pivot position in inertial frame
%    alpha         - pitch angle [rad]
%    x0dot, z0dot  - pivot velocity in inertial frame
%    alphadot      - pitch rate [rad/s]

% Position
x0 = -U_inf * t;
z0 = h_amp * sin(omega * t);

% Pitch angle
alpha = alpha0 + alpha_amp * sin(omega * t + phi_p);

% Velocities 
x0dot = -U_inf;
z0dot =  h_amp * omega * cos(omega * t);
alphadot =  alpha_amp * omega * cos(omega * t + phi_p);
end
