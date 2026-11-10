function dy = fun_ode_NL(T,X,MFF,N,Extr,Extr_tr,y,eps_Lp,eps_Lm,ssp,sfp,ssm,sfm,B,csi,Ea,Emp,Emm,EmT,b,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0)

dy = [MFF\(-fun_kff(X(2311:3060),X(3061:3810),X(61:810),X(811:1560),X(1561:2310),Extr_tr,Extr,N,csi,y,Ea,Emp,Emm,EmT,eps_Lp,eps_Lm,b,ssp,sfp,ssm,sfm,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,xisp0,xism0,xiT0)*(X(31:60)+X(1:30)/2000)) 
    X(1:30)
    (X(61:2310) - fun_xi(X(2311:3060),ssp,sfp,ssm,sfm,xisp0,xism0,xiT0))*10^3 % contiene tutti gli xi in colonna
    X(2311:3060) - (Ea + X(61:810)*(Emp-Ea) + X(811:1560)*(Emm-Ea) + X(1561:2310)*(EmT-Ea)) .* (X(3061:3810) - eps_Lp*X(61:810) - eps_Lm*X(811:1560)) %sigma
    (X(3061:3810) - fun_eps(B,N,X(31:60),Extr))*10^2 ];  % epsilon

end
