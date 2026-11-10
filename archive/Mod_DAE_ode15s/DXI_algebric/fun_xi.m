function xi = fun_xi(sigma,ssp,sfp,ssm,sfm,Dxisp,Dxism,DxiT)

global sigma_old xisp0 xism0 xiT0

indp = find(sigma>ssp & sigma < sfp & sigma > sigma_old);
indm = find(sigma<ssm & sigma > sfm & sigma < sigma_old);

xisp = xisp0; xism = xism0; xiT = xiT0;

if ~isempty(indp)
    xisp(indp) = (1-cos(pi*(sigma(indp)-ssp)/(sfp-ssp)))/2 + xisp0(indp) .* (1+cos(pi*(sigma(indp)-ssp)/(sfp-ssp)))/2;
    xism(indp) = xism0(indp) .* (1-xisp(indp))./(1-xisp0(indp));
    xiT(indp) = xiT0(indp) .* (1-xisp(indp))./(1-xisp0(indp));
    Dxisp(indp) = pi / (ssp-sfp) * sin(pi*(sigma(indp)-ssp)/(sfp-ssp)) .*(xisp0(indp)-1)/2;
    Dxism(indp) = xism0(indp) ./ (xisp0(indp) - 1) .* Dxisp(indp);
    DxiT(indp) = xiT0(indp) ./ (xisp0(indp) - 1) .* Dxisp(indp);
end
if ~isempty(indm)
    xism(indm) = (1-cos(pi*(sigma(indm)-ssm)/(sfm-ssm)))/2 + xism0(indm) .* (1+cos(pi*(sigma(indm)-ssm)/(sfm-ssm)))/2;
    xisp(indm) = xisp0(indm) .* (1-xism(indm))./(1-xism0(indm));
    xiT(indm) = xiT0(indm) .* (1-xism(indm))./(1-xism0(indm));
    Dxism(indm) = pi / (ssm-sfm) * sin(pi*(sigma(indm)-ssm)/(sfm-ssm)) .*(xism0(indm)-1)/2;
    Dxisp(indm) = xisp0(indm) ./ (xism0(indm) - 1) .* Dxism(indm);
    DxiT(indm) = xiT0(indm) ./ (xism0(indm) - 1) .* Dxism(indm);
end
xi = [xisp ; xism ; xiT ; Dxisp ; Dxism ; DxiT];

end
