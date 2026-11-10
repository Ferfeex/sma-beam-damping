function xdot = fun_ode(T,X,M,C,K,F,N)
xdot = [(-M\C)*X(1:3*N-4) - (M\K)*X(3*N-3:end) + M\(F*cos(5*T)) ; X(1:3*N-4)];
% xdot = [(-M\C)*X(1:3*N-4) - (M\K)*X(3*N-3:end) ; X(1:3*N-4)];
end
