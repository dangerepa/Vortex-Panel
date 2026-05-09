function [u, w] = vor2d(Gamma, x, z, xi, zi)
%VOR2D  Velocity induced by a 2-D point vortex.
% Aidan Bryce Sabin, May 2026
% From "Low-Speed Aerodynamics From Wing Theory to Panel Methods" by Katz
% and Plotkin
%  Computes the velocity (u, w) at field point (x, z) due to a vortex of
%  strength Gamma located at (xi, zi).
%
%  Based on the Biot-Savart law for a 2-D vortex filament:
%    V = Gamma / (2*pi*r),  directed perpendicular to r
%    u =  Gamma/(2*pi) *  (z-zi) / r^2
%    w = -Gamma/(2*pi) *  (x-xi) / r^2
%
%  A cut-off radius of 1e-6 is used to avoid singularity when the field
%  point coincides with the vortex location.
%
%  Inputs:
%    Gamma  - vortex strength [m^2/s]
%    x, z   - field point coordinates
%    xi, zi - vortex location coordinates
%
%  Outputs:
%    u  - x-component of induced velocity
%    w  - z-component of induced velocity

R_CUTOFF = 1e-6;

rx = x - xi;
rz = z - zi;
r = sqrt(rx^2 + rz^2);

if r < R_CUTOFF
    u = 0;
    w = 0;
    return
end

coeff = Gamma / (2*pi*r^2);
u = coeff * rz;
w = -coeff * rx;
end
