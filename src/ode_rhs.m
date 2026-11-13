function xdot = ode_rhs(T,X,Minv,N)
% Right-hand side of the ODE  M a + C v + K u = F  with the constant inverse mass matrix precomputed.
% State X = [velocities; displacements]; K, C and the load F are the globals KFF, CFF, FFF.
global KFF CFF FFF
n = 2*N-2;
v = X(1:n);  u = X(n+1:end);
if isempty(FFF), FFF = 0; end
xdot = [Minv*(-KFF*u - CFF*v + FFF); v];
end
