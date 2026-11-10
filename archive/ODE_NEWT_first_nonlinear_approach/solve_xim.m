function xism = solve_xim(sigma,xisp,ssp,sfp,ssm,sfm)
global sigma_old xisp0 xism0
xism = xism0;
indp = find(sigma>=sigma_old & sigma>ssp & sigma<sfp);
indm = find(sigma<=sigma_old & sigma<ssm & sigma>sfm);
if any(indm)
    xism(indm) = 1/2*(1-cos(pi*(sigma(indm)-ssm)/(sfm-ssm)))+(xism0(indm)/2).*(1+cos(pi*(sigma(indm)-ssm)/(sfm-ssm)));
end
if any(indp)
    xism(indp)=xism0(indp).*(1-xisp(indp))./(1-xisp0(indp));
end
end
