function xdot = fun_ode_1D(T,X,M,N,F)
global KFF CFF
xdot = [M\(-KFF*(X(2*N-1:end))-CFF*(X(1:2*N-2)))
    X(1:2*N-2)];
end
