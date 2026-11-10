function xi = fun_xi(sigma,ssp,sfp,ssm,sfm,xisp0,xism0,xiT0)
global sigma_old xisp_old xism_old xiT_old
% INDECES OF THOSE PARTS OF THE BEAM IN MD3+ AND MD3- PHASE DIAGRAM ZONES
indp = find(sigma>ssp & sigma < sfp & sigma >= sigma_old);
indm = find(sigma<ssm & sigma > sfm & sigma <= sigma_old);
% INITIALIZING XI TO X0
xisp = xisp_old; xism = xism_old; xiT = xiT_old;
% EVALUATION OF XI WHERE TRANSFORMATION IS TAKING PLACE
if any(indp)~=0
    xisp(indp) = (1-cos(pi*(sigma(indp)-ssp)/(sfp-ssp)))/2 + xisp0(indp).*(1+cos(pi*(sigma(indp)-ssp)/(sfp-ssp)))/2;
    xism(indp) = xism0(indp) .* (1-xisp(indp))./(1-xisp0(indp));
    xiT(indp) = xiT0(indp) .* (1-xisp(indp))./(1-xisp0(indp));
end
if any(indm)~=0
    xism(indm) = (1-cos(pi*(sigma(indm)-ssm)/(sfm-ssm)))/2 + xism0(indm) .* (1+cos(pi*(sigma(indm)-ssm)./(sfm-ssm)))/2;
    xisp(indm) = xisp0(indm) .* (1-xism(indm))./(1-xism0(indm));
    xiT(indm) = xiT0(indm) .* (1-xism(indm))./(1-xism0(indm));
end
xi = [xisp ; xism ; xiT];
end
