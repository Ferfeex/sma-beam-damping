function F = myfunsolve(X,epsilon,xisp0,xism0,xiT0,sigma_old,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm)
F = [ X(1:750) - (Ea + X(751:1500)*(Emp-Ea) + X(1501:2250)*(Emm-Ea) + X(2251:3000)*(EmT-Ea)) .* (epsilon - eps_Lp*X(751:1500) - eps_Lm*X(1501:2250))
      X(751:1500) - solve_xip(X(1:750),sigma_old,xisp0,xism0,X(1501:2250),ssp,sfp,ssm,sfm)
      X(1501:2250) - solve_xim(X(1:750),sigma_old,xisp0,xism0,X(751:1500),ssp,sfp,ssm,sfm)
      X(2251:3000) - solve_xit(X(1:750),X(751:1500),X(1501:2250),sigma_old,xiT0,xisp0,xism0,ssp,sfp,ssm,sfm) ];
end
