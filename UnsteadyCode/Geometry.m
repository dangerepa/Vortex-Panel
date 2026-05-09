function [cor, col_wf, vb_wf] = Geometry(c, N, I_t, cor, camber_x, camber_z)
%GEOMETRY  Define panel geometry in the wing-fixed frame.
% Aidan Bryce Sabin, May 2026
%  Works for both a flat plate (camber_x/z empty or not supplied) and a
%  general cambered airfoil (camber line supplied as x/z coordinates).
%
%  Panel layout :
%    Corner points   : N+1 points along the camber line
%    Bound vortex    : 1/4-chord of each panel (from its LE corner)
%    Collocation pt  : 3/4-chord of each panel
%    n_hat           : normal to the panel (outward from upper surface)
%    t_hat           : tangent along the panel (LE → TE direction)
%
%  Inputs:
%    c         - chord length [m]
%    N         - number of panels
%    I_t       - current time index
%    cor       - corner point structure (time history)
%    camber_x  - camber line x-coordinates in wing frame [1 x M], LE→TE
%    camber_z  - camber line z-coordinates in wing frame [1 x M]
%              If empty or omitted → flat plate (z = 0)
%
%  Outputs:
%    cor     - updated corner point structure (wing-fixed)
%    col_wf  - collocation points (wing-fixed): fields X, Z, nhat, that
%    vb_wf   - bound vortex positions (wing-fixed): fields X, Z

%  Flat plate fallback 
if nargin < 5 || isempty(camber_x)
    camber_x = linspace(0, c, N+1);
    camber_z = zeros(1, N+1);
end

%  Resample camber line to exactly N+1 evenly-spaced corner points 
% "Evenly spaced" here means equal arc-length intervals along the camber line.
arc = [0, cumsum(sqrt(diff(camber_x).^2 + diff(camber_z).^2))];
arc_total = arc(end);
arc_uniform = linspace(0, arc_total, N+1);

x_cor = interp1(arc, camber_x, arc_uniform, 'pchip');
z_cor = interp1(arc, camber_z, arc_uniform, 'pchip');

%  Store corner points 
cor(I_t).X = x_cor;
cor(I_t).Z = z_cor;

%  Panel tangent and normal vectors 
% Panel j runs from corner j to corner j+1.
% t_hat = unit vector from corner j to corner j+1
% n_hat = t_hat rotated 90 deg counter-clockwise (points to upper surface)
dx = diff(x_cor); % 1×N
dz = diff(z_cor); % 1×N
dl = sqrt(dx.^2 + dz.^2); % panel arc lengths

that = [dx; dz] ./ dl; % 2×N  unit tangent
nhat = [-dz; dx] ./ dl; % 2×N  unit normal 

%  Bound vortex positions: 1/4-chord of each panel 
x_vb = x_cor(1:N) + 0.25*dx;
z_vb = z_cor(1:N) + 0.25*dz;

vb_wf.X = x_vb;
vb_wf.Z = z_vb;

%  Collocation points: 3/4-chord of each panel 
x_col = x_cor(1:N) + 0.75*dx;
z_col = z_cor(1:N) + 0.75*dz;

col_wf.X = x_col;
col_wf.Z = z_col;
col_wf.nhat = nhat;
col_wf.that = that;
col_wf.dl = dl; 
end
