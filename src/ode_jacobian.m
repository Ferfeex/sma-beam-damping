function J = ode_jacobian(T,X,Minv,N)
% Analytic Jacobian of ode_rhs.m:  d xdot / d X = [ -Minv*C  -Minv*K ; I  0 ]
% (the system is linear for frozen K, C; K, C are the globals written by the OutputFcn).
global KFF CFF
n = 2*N-2;
J = [-Minv*CFF, -Minv*KFF; eye(n), zeros(n)];
end
