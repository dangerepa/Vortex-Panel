function A = Build_A(col_wf, vb_wf, vs, x0, z0, alpha, x_pivot, I_t, N)
%BUILD_A  Assemble the influence coefficient matrix.
% Aidan Bryce Sabin, May 2026
%  The system has N+1 unknowns: [Γ_1, ..., Γ_N, Γ_new]
%  Rows 1..N  : no-penetration BC at each collocation point (wing-fixed frame)
%  Row  N+1   : Kelvin's circulation theorem
%
%  Influence of bound vortex j on collocation point i is computed in the
%  WING-FIXED frame (velocities are evaluated there for the BC).
%
%  Inputs:
%    col_wf   - collocation points in wing-fixed frame
%    vb_wf    - bound vortex positions in wing-fixed frame
%    vs       - wake vortex structure (inertial), needed for new vortex position
%    x0,z0    - pivot position (inertial)
%    alpha    - pitch angle
%    x_pivot  - pivot chord location
%    I_t      - time index
%    N        - number of panels
%
%  Outputs:
%    A  - (N+1) × (N+1) influence matrix

A = zeros(N+1, N+1);

ca = cos(alpha);  sa = sin(alpha);
R_ItoW = [ca, sa; -sa, ca]; % Rotation: inertial -> wing-fixed

%  Rows 1:N: no-penetration BC 
for i = 1:N
    xi_wf = col_wf.X(i);
    zi_wf = col_wf.Z(i);
    ni = col_wf.nhat(:,i); % normal in wing-fixed frame

 % Influence of each bound vortex on collocation pt i
    for j = 1:N
        [u_wf, w_wf] = vor2d(1.0, xi_wf, zi_wf, vb_wf.X(j), vb_wf.Z(j));
        vel_wf = [u_wf; w_wf];
        A(i,j) = dot(vel_wf, ni);
    end

 % Influence of the new vortex on collocation pt i
 % New vortex position is in inertial frame -> transform to wing-fixed
    x_new_I = vs(I_t).X(end);
    z_new_I = vs(I_t).Z(end);

 % Transform new vortex to wing-fixed frame
    p_I = [x_new_I - x0; z_new_I - z0];
    p_rot = R_ItoW * p_I;
    x_new_wf = p_rot(1) + x_pivot;
    z_new_wf = p_rot(2);

    [u_wf, w_wf] = vor2d(1.0, xi_wf, zi_wf, x_new_wf, z_new_wf);
    vel_wf = [u_wf; w_wf];
    A(i, N+1) = dot(vel_wf, ni);
end

%  Row N+1: Kelvin
% Σ Γ_bound + Γ_new = const  ->  all coefficients = 1
A(N+1, :) = 1;
end
