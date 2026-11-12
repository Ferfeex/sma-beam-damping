function xdot = fun_ode_fast(T,X,Minv,N)
% Same right-hand side as fun_ode_1D.m but with the constant inverse mass matrix precomputed
% (no M\ solve at every call).  State X = [velocities; displacements].
global KFF CFF FFF
n = 2*N-2;
v = X(1:n);  u = X(n+1:end);
if isempty(FFF), FFF = 0; end
xdot = [Minv*(-KFF*u - CFF*v + FFF); v];
end
