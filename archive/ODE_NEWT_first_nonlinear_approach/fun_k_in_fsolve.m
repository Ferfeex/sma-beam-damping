function KFF = fun_k_in_fsolve(X,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm,b,N,Extr,Extr_tr,csi,y,B,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44)
global sigma_i xisp xism xiT epsilon sigma_old xisp0 xism0 xiT0 ind_old epsyylon xysp xysm
% DISPLACEMENT FIELD
XA = [X(31:60);zeros(2,1)];
X_k = zeros(4*(N-1),1);
for k = 1:N-1
    X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * XA;
end
epsilon = B * X_k;
% FIRST SIGMA GUESS
sigma_i = sigma_old + (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - epsyylon(:,end) - eps_Lp*(xisp-xysp(:,end)) - eps_Lm*(xism-xysm(:,end)));
% NEWTON RAPHSON + FSOLVE (WHEN NEWTON DOES NOT CONVERGE)
i = 0; % Iterations counter
index = 1; % Initialization of Newton's non-converged indeces
%%%%% LOCAL STATE ACCORDING TO MARTENSITIC PHASE DIAGRAM %%%%%%%%%%%%%%%%%%
indp = find(sigma_i > ssp & sigma_i < sfp & abs(sigma_i) > abs(sigma_old));
indp = setdiff(indp,ind_old);
indm = find(sigma_i < ssm & sigma_i > sfm & abs(sigma_i) > abs(sigma_old));
indm = setdiff(indm,ind_old);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dxisp = zeros(size(xisp)); Dxism = zeros(size(xism)); DxiT = zeros(size(xiT));
xisp = xisp0; xism = xism0; xiT = xiT0;
while ~isempty(index) && i < 500
    % MARTENSITIC FRACTIONS
    xisp(indp) = 1/2 * (1-cos(pi*(sigma_i(indp)-ssp)/(sfp-ssp))) + (xisp0(indp)/2) .* (1+cos(pi*(sigma_i(indp)-ssp)/(sfp-ssp)));
    xiT(indp) = xiT0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
    xism(indp) = xism0(indp) .* (1 - xisp(indp))./(1 - xisp0(indp));
    xism(indm) = 1/2 * (1-cos(pi*(sigma_i(indm)-ssm)/(sfm-ssm))) + (xism0(indm)/2) .* (1+cos(pi*(sigma_i(indm)-ssm)/(sfm-ssm)));
    xiT(indm) = xiT0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
    xisp(indm) = xisp0(indm) .* (1 - xism(indm))./(1 - xism0(indm));
    % SIGMA FOR COMPARISON
    sigmas = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - eps_Lp*xisp - eps_Lm*xism);
    % MARTENSITIC FRACTIONS DERIVATIVES
    Dxisp(indp) = 1/2*(xisp0(indp)-1) .* pi/(ssp-sfp) .* sin(pi*(sigma_i(indp)-ssp)/(sfp-ssp));
    DxiT(indp) = Dxisp(indp) .* (xiT0(indp)./(xisp0(indp)-1));
    Dxism(indp) = Dxisp(indp) .* (xism0(indp)./(xisp0(indp)-1));
    Dxism(indm) = 1/2*(xism0(indm)-1) .* pi/(ssm-sfm) .* sin(pi*(sigma_i(indm)-ssm)/(sfm-ssm));
    DxiT(indm) = Dxism(indm) .* (xiT0(indm)./(xism0(indm)-1));
    Dxisp(indm) = Dxism(indm) .* (xisp0(indm)./(xism0(indm)-1));
    % YOUNG'S MODULUS DERIVATIVE
    dE = Dxisp*(Emp-Ea) + Dxism*(Emm-Ea) + DxiT*(EmT-Ea);
    % NEWTON'S RESIDUAL, EVALUATION & DERIVATIVE
    F = sigma_i-sigmas;
    index = find(abs(F)>10^-6);
    dF = 1 - (dE.*(epsilon-eps_Lp*xisp-eps_Lm*xism) - (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)).*(-eps_Lp*Dxisp-eps_Lm*Dxism));
    % SIGMA UPDATE
    sigma_i(index) = sigma_i(index) - F(index)./dF(index);
    % UNION OF PHASE DIAGRAM & NON-CONVERGENCE INDECES
    indp = indp(ismember(indp,index));
    indm = indm(ismember(indm,index));
    i = i+1; % Counter update
end
% FSOLVE IN CASE OF NEWTON'S CONVERGENCE FAILURE
if i > 499
    % RE-INITIALIZATION OF THE VARIABLES WITH VALUES FROM NEWTON
    sigma_i = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon - eps_Lp*xisp - eps_Lm*xism);
    % INITIAL GUESS FOR FSOLVE
    X0 = [sigma_i ; xisp ; xism ; xiT];
    % JACOBIAN SPARSITY PATTERN FOR TRUST-REGION METHOD
    JB = zeros(3000);
    EY = repmat(eye(750),3,3);
    EY2 = repmat(eye(750),1,3);
    JB(1:2250,1:2250) = EY;
    JB(2251:3000,751:3000) = EY2;
    JB(1:750,2251:3000) = eye(750);
    options = optimoptions('fsolve','JacobPattern',JB,'Algorithm','trust-region');
    % NONLINEAR SPACE EQUATIONS SOLUTION 
    X = fsolve(@(X)myfunsolve_ode(X,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm) , X0, options);
    % SOLVED VARIABLES
    sigma_i = X(1:750);
    xisp = X(751:1500);
    xism = X(1501:2250);
    xiT = X(2251:3000);
    % PHASE DIAGRAM INDECES EVALUATION FOR XI DERIVATIVES
    indp = find(sigma_i > ssp & sigma_i < sfp & abs(sigma_i) > abs(sigma_old));
    indm = find(sigma_i < ssm & sigma_i > sfm & abs(sigma_i) > abs(sigma_old));
    % XI DERIVATIVES
    Dxisp = zeros(size(xisp)); Dxism = zeros(size(xism)); DxiT = zeros(size(xiT));
    Dxisp(indp) = 1/2*(xisp0(indp)-1) .* pi/(ssp-sfp) .* sin(pi*(sigma_i(indp)-ssp)/(sfp-ssp));
    DxiT(indp) = Dxisp(indp) .* (xiT0(indp)./(xisp0(indp)-1));
    Dxism(indp) = Dxisp(indp) .* (xism0(indp)./(xisp0(indp)-1));
    Dxism(indm) = 1/2*(xism0(indm)-1) .* pi/(ssm-sfm) .* sin(pi*(sigma_i(indm)-ssm)/(sfm-ssm));
    DxiT(indm) = Dxism(indm) .* (xiT0(indm)./(xism0(indm)-1));
    Dxisp(indm) = Dxism(indm) .* (xisp0(indm)./(xism0(indm)-1));
end
% NEW YOUNG'S MODULUS & DERIVATIVE
E = Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea);
dE = Dxisp*(Emp-Ea) + Dxism*(Emm-Ea) + DxiT*(EmT-Ea);
% TANGENT MATRIX & STIFFNESS MATRIX FOR K-TH ELEMENT
C_tot = E ./ ( 1-(epsilon-eps_Lp*xisp-eps_Lm*xism).*dE + E.*(eps_Lp*Dxisp+eps_Lm*Dxism) );
Ck = reshape(C_tot,length(y),length(csi),N-1);
Kk = zeros(4,4,N-1);
for k = 1:N-1
    Kk(:,:,k) = b*[trapz(csi,trapz(-y,(Ck(:,:,k).*BB11))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB12))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB13))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB14)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB21))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB22))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB23))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB24)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB31))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB32))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB33))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB34)))
        trapz(csi,trapz(-y,(Ck(:,:,k).*BB41))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB42))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB43))) , trapz(csi,trapz(-y,(Ck(:,:,k).*BB44)))];
end
% STIFFNESS MATRIX ASSEMBLY
K = zeros(2*N,2*N,N-1);
for k = 1:N-1
    K(:,:,k) = Extr_tr(:,:,k) * Kk(:,:,k) * Extr(:,:,k);
end
K = sum(K,3);
KFF = K(1:2*N-2,1:2*N-2);
end
