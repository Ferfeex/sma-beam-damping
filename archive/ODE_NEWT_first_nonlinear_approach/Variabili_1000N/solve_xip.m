function xisp = solve_xip(sigma,sigma_old,xisp0,xism0,xism,ssp,sfp,ssm,sfm)
xisp = xisp0;
indp = find(sigma>=sigma_old & sigma>ssp & sigma<sfp);
indm = find(sigma<=sigma_old & sigma<ssm & sigma>sfm);
if any(indp)
    xisp(indp) = 1/2*(1-cos(pi*(sigma(indp)-ssp)/(sfp-ssp)))+(xisp0(indp)/2).*(1+cos(pi*(sigma(indp)-ssp)/(sfp-ssp)));
end
if any(indm)
    xisp(indm)=xisp0(indm).*(1-xism(indm))./(1-xism0(indm));
end
end
