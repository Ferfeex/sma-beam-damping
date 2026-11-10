function xiT = solve_xit(sigma,xisp,xism,sigma_old,xiT0,xisp0,xism0,ssp,sfp,ssm,sfm)
xiT = xiT0;
indp = find(sigma>=sigma_old & sigma>ssp & sigma<sfp);
indm = find(sigma<=sigma_old & sigma<ssm & sigma>sfm);
if any(indp)
    xiT(indp) = xiT0(indp).*(1-xisp(indp))./(1-xisp0(indp));
end
if any(indm)
    xiT(indm) = xiT0(indm).*(1-xism(indm))./(1-xism0(indm));
end
end
