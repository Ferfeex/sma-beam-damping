function xdot = fun_ode_1D(T,X,M,F,N,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm,b,Extr,Extr_tr,csi,y,B,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44)
xdot = [M\(-fun_k_in_fsolve(X,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm,b,N,Extr,Extr_tr,csi,y,B,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44)*(X(2*N-1:end)+X(1:2*N-2)/2300))
    X(1:2*N-2)];
clc
fprintf('T = %s \n',T)
end
