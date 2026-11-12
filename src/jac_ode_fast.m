function J = jac_ode_fast(T,X,Minv,N)
% Analytic Jacobian of fun_ode_1D / fun_ode_fast:  d xdot / d X = [ -Minv*C  -Minv*K ; I  0 ]
% (the system is linear for frozen K, C; K, C are the globals written by the OutputFcn).
global KFF CFF
n = 2*N-2;
J = [-Minv*CFF, -Minv*KFF; eye(n), zeros(n)];
end
