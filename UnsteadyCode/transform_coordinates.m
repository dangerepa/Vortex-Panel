function [col, vb] = transform_coordinates(col_wf, vb_wf, x0, z0, alpha, ...
                                             x_pivot, I_t, col, vb)
%TRANSFORM_COORDINATES  Map wing-fixed frame -> inertial frame.
% Aidan Bryce Sabin, May 2026
%  The wing-fixed frame is rotated by alpha and translated so that the
%  pivot point maps to (x0, z0) in the inertial frame.
%
%  Rotation matrix  R_{W->I}:
%    [x_I]   [cos(alpha) -sin(alpha)] [x_w - x_pivot]   [x0]
%    [z_I] = [sin(alpha) cos(alpha)] [z_w] + [z0]
%
%  Note: x_pivot shifts the chord origin so that rotation is about the
%  correct pivot location along the chord.
%
%  Inputs:
%    col_wf   - collocation points in wing-fixed frame
%    vb_wf    - bound vortex positions in wing-fixed frame
%    x0, z0   - pivot position in inertial frame
%    alpha    - pitch angle [rad]
%    x_pivot  - pivot x-location along chord (wing frame)
%    I_t      - current time index
%    col, vb  - structures to update
%
%  Outputs:
%    col  - updated collocation point structure (inertial)
%    vb   - updated bound vortex structure (inertial)

ca = cos(alpha);
sa = sin(alpha);
R = [ca, -sa; sa, ca]; % Rotation: wing -> inertial

%  Collocation points 
N = length(col_wf.X);
col_wf_mat = [col_wf.X - x_pivot; col_wf.Z]; % shift pivot to origin
col_I = R * col_wf_mat + [x0; z0]; % rotate + translate

col(I_t).X = col_I(1,:);
col(I_t).Z = col_I(2,:);

% Transform normal and tangent vectors 
col(I_t).nhat = R * col_wf.nhat; % 2×N in inertial frame
col(I_t).that = R * col_wf.that;

%  Bound vortex positions 
vb_wf_mat = [vb_wf.X - x_pivot; vb_wf.Z];
vb_I = R * vb_wf_mat + [x0; z0];

vb(I_t).X = vb_I(1,:);
vb(I_t).Z = vb_I(2,:);
% Gamma field initialised to zero; filled after solve
if isempty(vb(I_t).Gamma)
    vb(I_t).Gamma = zeros(1, N);
end
end
