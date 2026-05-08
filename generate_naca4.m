function [x, y] = generate_naca4(naca_str, N, filename)
% GENERATE_NACA4  Generate a NACA 4-digit airfoil profile
%
%   [x, y] = generate_naca4(naca_str, N, filename)
%   Inputs:
%     naca_str - 4-digit NACA string, e.g. '2412' or '0018'
%     N - points per surface, excluding trailing edge (default: 100)
%     filename  - (optional) Selig .dat output path
%
%   Outputs:
%     x, y - coordinates (TE -> upper -> LE -> lower -> TE)
%
%   Examples:
%     [x, y] = generate_naca4('2412', 100);
%     [x, y] = generate_naca4('0018', 100, 'NACA0018.dat');
% This uses the formula from https://www.pdas.com/naca456thick4.html
    if nargin < 2 || isempty(N), N = 100; end
    if nargin < 3, filename = ''; end

    % 4-digit designation 
    assert(numel(naca_str) == 4, 'Input must be a 4-character string, e.g. ''2412''.');
    digits = naca_str - '0'; % char to int array
    m = digits(1) / 100; % max camber
    p = digits(2) / 10; % location of max camber
    t = (digits(3)*10 + digits(4)) / 100;  % max thickness

    % Cosine-clustered x stations: TE -> LE
    beta  = linspace(0, pi, N+1);
    x_cos = 0.5 * (1 - cos(beta)); % 0 -> 1, LE to TE

    % Thickness distribution
    yt = 5*t .* (  0.2969 * sqrt(x_cos) ...
                 - 0.1260 * x_cos ...
                 - 0.3516 * x_cos.^2 ...
                 + 0.2843 * x_cos.^3 ...
                 - 0.1015 * x_cos.^4);

    % Camber line and its gradient 
    yc = zeros(size(x_cos));
    dyc = zeros(size(x_cos));

    if m > 0 && p > 0
        % Forward of max camber point
        idx1 = x_cos <= p;
        yc(idx1) = (m / p^2) .* (2*p*x_cos(idx1) - x_cos(idx1).^2);
        dyc(idx1) = (2*m / p^2) .* (p - x_cos(idx1));

        % Aft of max camber point
        idx2 = ~idx1;
        yc(idx2) = (m / (1-p)^2) .* (1 - 2*p + 2*p*x_cos(idx2) - x_cos(idx2).^2);
        dyc(idx2) = (2*m / (1-p)^2) .* (p - x_cos(idx2));
    end
    % (symmetric: yc and dyc stay zero)

    theta = atan(dyc); % surface normal angle

    % Upper and lower surface coordinates 
    xu = x_cos - yt .* sin(theta);
    yu = yc + yt .* cos(theta);

    xl = x_cos + yt .* sin(theta);
    yl = yc - yt .* cos(theta);

    % Format: TE -> upper -> LE -> lower -> TE
    x = [fliplr(xu), xl(2:end)]';
    y = [fliplr(yu), yl(2:end)]';

    % Optional: write .dat file 
    if ~isempty(filename)
        header = sprintf('NACA %s Airfoil', naca_str);
        fid = fopen(filename, 'w');
        if fid == -1, error('Could not open file: %s', filename); end
        fprintf(fid, '%s\n', header);
        fprintf(fid, '  %.6f  %.6f\n', [x, y]');
        fclose(fid);
        fprintf('Written to %s  (%d points)\n', filename, length(x));
    end
end