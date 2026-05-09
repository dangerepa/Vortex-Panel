function [camber_x, camber_z] = load_airfoil(source, c, varargin)
%LOAD_AIRFOIL  Extract camber line from an airfoil definition.
% Aidan Bryce Sabin, May 2026
%  Three usage modes:
%
%  1) NACA 4-digit string (generated analytically):
%       [cx, cz] = load_airfoil('NACA2412', c)
%       [cx, cz] = load_airfoil('NACA0012', c)
%
%  2) .dat file path (Selig or Lednicer format, upper+lower surface):
%       [cx, cz] = load_airfoil('path/to/airfoil.dat', c)
%
%  3) Raw coordinate matrix [x, z] already loaded (upper+lower surface rows):
%       [cx, cz] = load_airfoil(coords_matrix, c)
%
%  In modes 2 and 3 the camber line is extracted by averaging upper and
%  lower surface z-values at the same x stations.
%
%  All coordinates are scaled to chord length c and normalised so that
%  the leading edge is at x=0.
%
%  Inputs:
%    source  - NACA string, .dat filepath, or Nx2 coordinate matrix
%    c - chord length to scale to [m]
%    varargin(1) - number of points along camber line (default 200)
%
%  Outputs:
%    camber_x - camber line x-coordinates [1×M], LE->TE, range [0, c]
%    camber_z - camber line z-coordinates [1×M]

M = 200;
if ~isempty(varargin), M = varargin{1}; end

%  Mode 1: NACA 4-digit 
if ischar(source) && strncmpi(source, 'NACA', 4)
    digits = source(5:end);
    if length(digits) ~= 4
        error('load_airfoil: only 4-digit NACA supported, e.g. NACA2412');
    end
    m = str2double(digits(1)) / 100; % max camber
    p = str2double(digits(2)) / 10; % location of max camber
 % thickness digit unused for camber line

    x_c = linspace(0, 1, M); % normalised x
    z_c = zeros(1, M);

    if m > 0 && p > 0
        for i = 1:M
            xn = x_c(i);
            if xn < p
                z_c(i) = (m/p^2) * (2*p*xn - xn^2);
            else
                z_c(i) = (m/(1-p)^2) * (1 - 2*p + 2*p*xn - xn^2);
            end
        end
    end

    camber_x = x_c * c;
    camber_z = z_c * c;
    fprintf('load_airfoil: NACA %s camber line generated (%d points)\n', digits, M);
    return
end

%  Mode 2: .dat file 
if ischar(source)
    raw = load_dat_file(source);
else
    raw = source; % Mode 3: already a matrix
end

% Normalise x to [0,1] range
raw(:,1) = (raw(:,1) - min(raw(:,1))) / (max(raw(:,1)) - min(raw(:,1)));

%  Split upper and lower surfaces 
%  starts at TE (x=1), goes over upper surface to LE (x=0),
% then back along lower to TE. Find the LE index (minimum x).
[~, ile] = min(raw(:,1));

upper = raw(1:ile, :); % TE -> LE (upper)
lower = raw(ile:end, :); % LE -> TE (lower)

% Ensure both go LE -> TE for interpolation
if upper(1,1) > upper(end,1),  upper = flipud(upper); end
if lower(1,1) > lower(end,1),  lower = flipud(lower); end

%  Camber line:
x_common = linspace(0, 1, M);
z_upper = interp1(upper(:,1), upper(:,2), x_common, 'pchip', 'extrap');
z_lower = interp1(lower(:,1), lower(:,2), x_common, 'pchip', 'extrap');
z_camber = 0.5 * (z_upper + z_lower);

camber_x = x_common * c;
camber_z = z_camber * c;
fprintf('load_airfoil: camber line extracted from coordinates (%d points)\n', M);
end

%
function coords = load_dat_file(filepath)
    fid = fopen(filepath, 'r');
    if fid < 0
        error('load_airfoil: cannot open file %s', filepath);
    end

    coords = [];
    header_done = false;

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || line(1) == '#', continue; end

        nums = sscanf(line, '%f %f');
        if length(nums) == 2
       % Skip header line
            if ~header_done && (nums(1) > 10 || nums(2) > 10)
                header_done = true;
                continue
            end
            header_done = true;
            coords(end+1, :) = nums'; 
        end
    end
    fclose(fid);

    if isempty(coords)
        error('load_airfoil: no numeric coordinate pairs found in %s', filepath);
    end
end
