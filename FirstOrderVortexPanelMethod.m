%First Order Vortex Panel Method, Aidan Bryce Sabin, May 2026
clc; clear; close all;

%% Inputs and stuff
Vinf = 1.0;
alpha_deg = 5;
alpha = deg2rad(alpha_deg);
%% Generate airfoil and assign variables
[x, y] = generate_naca4('2412', 100);
x1 = x(1:end-1);
y1 = y(1:end-1);
x2 = x(2:end);
y2 = y(2:end);

dx = x2 - x1;
dy = y2 - y1;

S = sqrt(dx.^2 + dy.^2);
phi = atan2(dy, dx);
xc = 0.5*(x1 + x2);
yc = 0.5*(y1 + y2);

n = length(S);

%% solve 2
warning('off', 'all');
A = zeros(n,n);
b = zeros(n,1);

nx = sin(phi);
ny = -cos(phi);
% Method from "Anderson's Fundamentals of Aerodynamics", 7th ed ch. 4.10
for i = 1:n-1 %if broken go back to n
    b(i) = -Vinf*(cos(alpha)*nx(i) + sin(alpha)*ny(i));

    for j = 1:n
        if i==j
            A(i,j)=0;
        else
            integrand = @(s) dtheta_dn_integrand( ...
                xc(i), yc(i), nx(i), ny(i), ...
                x1(j), y1(j), phi(j), s);

            A(i,j) = -integral(integrand, 0, S(j), ...
                'AbsTol',1e-10,'RelTol',1e-8) / (2*pi);
        end
    end
end

A(n,:) = 0;
A(n,1) = 1;
A(n,n) = 1;
b(n) = 0;

%solve
gamma = A\b;
gamma_avg = 0.5*(gamma(1:end-1) + gamma(2:end));
% Circulation
Gamma = sum(gamma .* S);

fprintf('alpha = %.2f deg\n', alpha_deg);
fprintf('Gamma = %.6f\n', Gamma);

%% Cp and velocity
Vt = gamma; % velocity just outside surface
Cp = 1 - (Vt./Vinf).^2;

Vt_avg = gamma_avg;
xc_avg = 0.5*(xc(1:end-1) + xc(2:end));
Cp_avg = 1 - (gamma_avg./Vinf).^2;

x_le = min(x);
c = max(x) - min(x);

% Reference points for moment
x_ref_LE = x_le;
y_ref_LE = 0;

x_ref_c4 = x_le + 0.25*c;
y_ref_c4 = 0;

% Differential force coefficients on each panel
dCx = -Cp .* nx .* (S/c);
dCy = -Cp .* ny .* (S/c);

% Differential moment coefficients about LE and quarter-chord
dCm_LE = ((xc - x_ref_LE).*dCy - (yc - y_ref_LE).*dCx)/c;
dCm_c4 = ((xc - x_ref_c4).*dCy - (yc - y_ref_c4).*dCx)/c;

% Total moment coefficients
Cm_LE = sum(dCm_LE);
Cm_c4 = sum(dCm_c4);

%These two should be close
CL_kj = 2*Gamma/(Vinf*c); %kutta joukowski lift
CL_p = sum(-Cp .* ny .* (S/c)); %pressure lift

fprintf('CL from Kutta-Joukowski = %.6f\n', CL_kj);
fprintf('CL from pressure = %.6f\n', CL_p);

fprintf('Cm_LE  = %.6f\n', Cm_LE);
fprintf('Cm_c/4 = %.6f\n', Cm_c4);


% Plot
figure; %Cp
yyaxis left
plot(xc_avg, Cp_avg, 'bo-');
ylabel('C_p');
set(gca,'YDir','reverse');
ylim([-2 1])

yyaxis right
plot(x, y, 'k-');
ylim([-0.4 0.4]);
ylabel('y');

xlabel('x');
title(sprintf('C_p Distribution and Airfoil Shape, \\alpha = %.1f^\\circ', alpha_deg));
grid on;
legend('C_p','Airfoil','Location','best');

figure; %control points
plot(x, y, 'k-'); hold on;
plot(xc, yc, 'ro');
axis equal;
xlabel('x'); ylabel('y');
title('Panel Geometry and Control Points');
grid on;

figure; %gamma
plot(1:n, gamma, 'o-');
xlabel('Panel index');
ylabel('\gamma');
title('Constant vortex strength on each panel');
grid on;

function val = dtheta_dn_integrand(xc, yc, nx_i, ny_i, xj1, yj1, phi_j, s)
xj = xj1 + s*cos(phi_j);
yj = yj1 + s*sin(phi_j);

dx = xc - xj;
dy = yc - yj;

val = (-dy*nx_i + dx*ny_i) ./ (dx.^2 + dy.^2);
end
