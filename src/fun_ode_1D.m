function xdot = fun_ode_1D(T,X,M,N)
% Right-hand side of the first-order system  X = [velocities; displacements]:
%   M*a = -K*u - C*v + F      (K, C change at run time, written by myfun.m)
%   u'  = v
% The load term FFF is zero except in the "pulse" cases. (The original 4-argument
% version had no force term; the force term FFF comes from Loss/Transient/fun_ode_1D.m.)
global KFF CFF FFF
if isempty(FFF), FFF = 0; end
xdot = [M\(-KFF*(X(2*N-1:end))-CFF*(X(1:2*N-2))+FFF)
    X(1:2*N-2)];
end
