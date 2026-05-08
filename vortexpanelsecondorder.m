%Second Order Vortex Panel Method, Aidan Sabin, May 2026
clc; clear; close all;
%% Inputs
[x, y] = generate_naca4('2412', 100); %100 makes 200 panels which is a good number
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

n = length(S); % number of panels
nn = n + 1; % number of nodal gamma unknowns

nx = sin(phi);
ny = -cos(phi);
Vinf = 1.0;
alpha_deg = 15;
alpha = deg2rad(alpha_deg);
%% Assemble n collocation equations for n+1 (nn) unknowns
A = zeros(nn, nn);
b = zeros(nn, 1);

for i = 1:n
    b(i) = -Vinf*(cos(alpha)*nx(i) + sin(alpha)*ny(i));

    for j = 1:n
        if i == j
            % principal value self-term for normal velocity on a straight
            % vortex panel is zero
            aL = 0;
            aR = 0;
        else
            integrandL = @(s) (1 - s./S(j)) .* dtheta_dn_integrand( ...
                xc(i), yc(i), nx(i), ny(i), x1(j), y1(j), phi(j), s);

            integrandR = @(s) (s./S(j)) .* dtheta_dn_integrand( ...
                xc(i), yc(i), nx(i), ny(i), x1(j), y1(j), phi(j), s);

            aL = -integral(integrandL, 0, S(j), ...
                'AbsTol',1e-10,'RelTol',1e-8) / (2*pi);

            aR = -integral(integrandR, 0, S(j), ...
                'AbsTol',1e-10,'RelTol',1e-8) / (2*pi);
        end

        % panel j uses nodal unknowns j and j+1
        A(i,j)   = A(i,j)   + aL;
        A(i,j+1) = A(i,j+1) + aR;
    end
end

%% Kutta 
% Upper and lower TE nodal strengths sum to zero
A(nn,:) = 0;
A(nn,1) = 1;
A(nn,nn) = 1; %sometimes this needs to be -1
b(nn) = 0;

gamma_node = A\b; %solve system

%% Panel-averaged gamma for plotting / circulation
gamma_panel = 0.5*(gamma_node(1:end-1) + gamma_node(2:end));

Gamma = sum(gamma_panel .* S);
% 
c = max(x) - min(x);
% CL = 2*Gamma/(Vinf*c);
%Gamma = sum(gamma_panel(2:end-1) .* S(2:end-1));
CL = 2*Gamma/(Vinf*c);
%print results
fprintf('alpha = %.2f deg\n', alpha_deg);
fprintf('Gamma = %.6f\n', Gamma);
fprintf('CL = %.6f\n', CL);

%% Pressure and Moment dist.
tx = cos(phi);
ty = sin(phi);

gamma_cp = 0.5*(gamma_node(1:end-1) + gamma_node(2:end));

Vt = zeros(n,1);

for i = 1:n
    % freestream tangential component
    Vt(i) = Vinf*(cos(alpha)*tx(i) + sin(alpha)*ty(i));

    for j = 1:n
        if i == j
            bL = 0;
            bR = 0;
        else
            integrandL = @(s) (1 - s./S(j)) .* dtheta_dt_integrand( ...
                xc(i), yc(i), tx(i), ty(i), x1(j), y1(j), phi(j), s);

            integrandR = @(s) (s./S(j)) .* dtheta_dt_integrand( ...
                xc(i), yc(i), tx(i), ty(i), x1(j), y1(j), phi(j), s);

            bL = -integral(integrandL, 0, S(j), ...
                'AbsTol',1e-10,'RelTol',1e-8)/(2*pi);

            bR = -integral(integrandR, 0, S(j), ...
                'AbsTol',1e-10,'RelTol',1e-8)/(2*pi);
        end

        Vt(i) = Vt(i) + bL*gamma_node(j) + bR*gamma_node(j+1);
    end

    % jump term: velocity just outside the surface
    Vt(i) = Vt(i) - 0.5*gamma_cp(i);
end


%% plot Cp
Cp = 1 - (Vt./Vinf).^2;
yyaxis left 
plot(xc, Cp, 'bo-'); 
ylabel('C_p'); 
set(gca,'YDir','reverse'); 
ylim([-2 1])

yyaxis right
plot(x, y, 'k-','LineWidth',1.2);
ylabel('y');
ylim([-0.4 0.4]);

xlabel('x');
title(sprintf('C_p Distribution and Airfoil Shape, \\alpha = %.1f^\\circ', alpha_deg));
grid on;
legend('C_p','Airfoil','Location','best');

%% Plot nodal gamma
figure;
plot(1:nn, gamma_node, 'o-');
xlim([0 201])
xlabel('Boundary-point index');
ylabel('\gamma node');
title('Nodal vortex-sheet strength (linear on each panel)');
grid on;

%% Plot panel-average gamma
figure;
plot(1:n, gamma_panel, 'o-');
xlabel('Panel index');
ylabel('\gamma panel avg');
title('Panel-averaged \gamma');
grid on;
%% Moment
x_le = min(x);
c = max(x) - min(x);

x_ref_LE = x_le;
y_ref_LE = 0;

x_ref_c4 = x_le + 0.25*c;
y_ref_c4 = 0;

dCx = -Cp .* nx .* (S/c);
dCy = -Cp .* ny .* (S/c);

dCm_LE = ((xc - x_ref_LE).*dCy - (yc - y_ref_LE).*dCx)/c;
dCm_c4 = ((xc - x_ref_c4).*dCy - (yc - y_ref_c4).*dCx)/c;

Cm_LE = sum(dCm_LE);
Cm_c4 = sum(dCm_c4);

fprintf('Cm_LE  = %.6f\n', Cm_LE);
fprintf('Cm_c/4 = %.6f\n', Cm_c4);
%% Cl
Gamma = sum(gamma_cp .* S);
CL_kj = 2*Gamma/(Vinf*c);

CL_p = sum(-Cp .* ny .* (S/c));

fprintf('CL from Kutta-Joukowski = %.6f\n', CL_kj);
fprintf('CL from pressure = %.6f\n', CL_p);
%% Functions
function val = dtheta_dn_integrand(xc, yc, nx_i, ny_i, xj1, yj1, phi_j, s)
    xj = xj1 + s*cos(phi_j);
    yj = yj1 + s*sin(phi_j);

    dx = xc - xj;
    dy = yc - yj;

    val = (-dy*nx_i + dx*ny_i) ./ (dx.^2 + dy.^2);
end
function val = dtheta_dt_integrand(xc, yc, tx_i, ty_i, xj1, yj1, phi_j, s)
    xj = xj1 + s*cos(phi_j);
    yj = yj1 + s*sin(phi_j);

    dx = xc - xj;
    dy = yc - yj;

    val = (-dy*tx_i + dx*ty_i) ./ (dx.^2 + dy.^2);
end