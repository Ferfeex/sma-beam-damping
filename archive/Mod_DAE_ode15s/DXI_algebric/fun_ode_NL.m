function ydot = fun_ode_NL(T,X,MFF,CFF,N,F,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,Lk,b)
% Mass matrix ---> (3810 x 3810)
% Valid for N=16, length(y)=10 and length(xe)=5

ydot = [ MFF\(-CFF*X(1:30) - fun_kff(X(5311:6060),X(61:810),X(811:1560),X(1561:2310),Extr_tr,Extr,N,csi,y,Lk,Ea,Emp,Emm,EmT,eps_Lp,eps_Lm,b,X(2311:3060),X(3061:3810),X(3811:4560))*X(31:60) + 0*F*cos(T)) %XF_DDOT
    X(1:30) %XF_DOT
    X(61:4560) - fun_xi(X(4561:5310),ssp,sfp,ssm,sfm,X(2311:3060),X(3061:3810),X(3811:4560)) %XI
    X(4561:5310) - (Ea + X(61:810) * (Emp - Ea) + X(811:1560) * (Emm - Ea) + X(1561:2310) * (EmT - Ea)) .* (X(5311:6060) - eps_Lp * X(61:810) - eps_Lm * X(811:1560)) %SIGMA
    X(5311:6060) - fun_eps(B,N,X(31:60),Extr) ];  %EPSILON
end
