function Gamma_vec = Solve_system(A, RHS)
%SOLVE_SYSTEM  Solve the linear system [A]{Γ} = {RHS}.
% Aidan Bryce Sabin, May 2026
%
%  Inputs:
%    A    - (N+1) × (N+1) influence matrix
%    RHS  - (N+1) × 1 right-hand side vector
%
%  Outputs:
%    Gamma_vec  - (N+1) × 1 vector of vortex strengths
%                 Gamma_vec(1:N)  = bound vortex strengths Γ_1..Γ_N
%                 Gamma_vec(N+1)  = newly shed vortex strength Γ_new

Gamma_vec = A \ RHS;
end
