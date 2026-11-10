function xi = functxi(sigma,ssp,sfp,ssm,sfm)
indp = find(sigma>ssp & sigma < sfp);
indm = find(sigma<ssm & sigma > sfm);

xi_sp = zeros(size(sigma)); xi_sm = zeros(size(sigma)); xi_T = ones(size(sigma));

if ~isempty(indp)
    xi_sp(indp) = (1-cos(pi*(sigma(indp)-ssp)/(sfp-ssp)))/2;
    xi_T(indp) = xi_T(indp) .* 1-xi_sp(indp);
end
if ~isempty(indm)
    xi_sm(indm) = (1-cos(pi*(sigma(indm)-ssm)/(sfm-ssm)))/2;
    xi_T(indm) = xi_T(indm) .* 1-xi_sm(indm);
end
xi = [xi_sp ; xi_sm ; xi_T];

end
