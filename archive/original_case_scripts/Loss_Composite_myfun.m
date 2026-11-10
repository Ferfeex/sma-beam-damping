function status = myfun(T,X,flag,Ea,Emp,Emm,EmT,ssp,sfp,ssm,sfm,eps_Lp,eps_Lm,b,Extr,Extr_tr,csi,y,B,N,BB11,BB12,BB13,BB14,BB21,BB22,BB23,BB24,BB31,BB32,BB33,BB34,BB41,BB42,BB43,BB44,p,KFF_HLE,CFF_HLE)
global KFF CFF tt yy xisp xism xiT xisp0 xism0 xiT0 indp_old indm_old sigma_old s p1 m1
global epsyylon xysp xysm xyT LF
if nargin < 3 || isempty(flag)
    % VARIABLES UPDATE AND EVALUATION
    tt = [tt T];
    yy = [yy X];
    % DISPLACEMENT FIELD
    XA = [X(31:60);zeros(2,1)];
    X_k = zeros(4*(N-1),1);
    for k = 1:N-1
        X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * XA;
    end
    epsilon = B * X_k;
    % FIRST SIGMA GUESS
    sigma_i = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon-eps_Lp*xisp-eps_Lm*xism);
    % NEWTON RAPHSON + FSOLVE (WHEN NEWTON DOES NOT CONVERGE)
    i = 0; % Iterations counter
    %%%%% LOCAL STATE ACCORDING TO MARTENSITIC PHASE DIAGRAM %%%%%%%%%%%%%%%%%%
    indp = find(sigma_i > ssp & sigma_i < sfp & sigma_i >= sigma_old);
    indm = find(sigma_i < ssm & sigma_i > sfm & sigma_i <= sigma_old);
    %%%%%%%%%%%%%%%%%%% UPDATING XI_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ex_indp = setdiff(indp_old, indp); ex_indm = setdiff(indm_old, indm);
    
    am = find(sigma_i > 0);
    ap = find(sigma_i < 0);
    if any(am)~=0
        m1 = setdiff(m1,am);
    end
    if any(ap)~=0
        p1 = setdiff(p1,ap);
    end
    p1 = union(p1,ex_indp);
    m1 = union(m1, ex_indm);
    indp = setdiff(indp,ex_indp); indp = setdiff(indp,p1);
    indm = setdiff(indm,ex_indm); indm = setdiff(indm,m1);
    
    indp_old = indp; indm_old = indm;
    
    sigmas = sigma_i;
    dF = zeros(750,1);
    Dxisp = zeros(size(xisp)); Dxism = zeros(size(xism)); DxiT = zeros(size(xiT));
    while ~isempty(indp) && ~isempty(indm)
        % MARTENSITIC FRACTIONS
        xisp(indp) = 1/2 * (1-cos(pi*(sigma_i(indp)-ssp)/(sfp-ssp))) + (xisp0(indp)/2) .* (1+cos(pi*(sigma_i(indp)-ssp)/(sfp-ssp)));
        xiT(indp) = xiT0(indp) .* ((1 - xisp(indp))./(1 - xisp0(indp)));
        xism(indp) = xism0(indp) .* ((1 - xisp(indp))./(1 - xisp0(indp)));
        xism(indm) = 1/2 * (1-cos(pi*(sigma_i(indm)-ssm)/(sfm-ssm))) + (xism0(indm)/2) .* (1+cos(pi*(sigma_i(indm)-ssm)/(sfm-ssm)));
        xiT(indm) = xiT0(indm) .* ((1 - xism(indm))./(1 - xism0(indm)));
        xisp(indm) = xisp0(indm) .* ((1 - xism(indm))./(1 - xism0(indm)));
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
        F = sigma_i-sigmas; % At those indeces ~=indp,indm this difference =0
        index = find(abs(F)>10^-6);
        dF(index) = 1 - (dE(index).*(epsilon(index)-eps_Lp*xisp(index)-eps_Lm*xism(index)) + (Ea + xisp(index)*(Emp-Ea) + xism(index)*(Emm-Ea) + xiT(index)*(EmT-Ea)).*(-eps_Lp*Dxisp(index)-eps_Lm*Dxism(index)));
        % SIGMA UPDATE
        sigma_i(index) = sigma_i(index) - F(index)./dF(index);
        % UNION OF PHASE DIAGRAM & NON-CONVERGENCE INDECES
        indp = indp(ismember(indp,index));
        indm = indm(ismember(indm,index));
        i = i+1; % Counter update
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST NEWTON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigma_i = sigmas;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    KFF = K(1:2*N-2,1:2*N-2) + KFF_HLE;
    epsyylon = [epsyylon, epsilon];
    xysp = [xysp , xisp];
    xysm = [xysm , xism];
    xyT = [xyT, xiT];
    s = [s , sigma_i];
    sigma_old = sigma_i;
    fprintf('t = %f [s] ', T)
    fprintf('i = %i [s] \n', i)
    if any(ex_indp)~=0 || any(ex_indm)~=0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW DAMPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LF(ex_indp) = polyval(p,epsilon(ex_indp));
        LF(ex_indm) = polyval(p,abs(epsilon(ex_indm)));
        C_d = C_tot .* LF / (2*pi*14.5);
        C_dk = reshape(C_d,length(y),length(csi),N-1);
        C_k = zeros(4,4,N-1);
        for k = 1:N-1
            C_k(:,:,k) = b*[trapz(csi,trapz(-y,(C_dk(:,:,k).*BB11))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB12))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB13))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB14)))
                trapz(csi,trapz(-y,(C_dk(:,:,k).*BB21))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB22))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB23))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB24)))
                trapz(csi,trapz(-y,(C_dk(:,:,k).*BB31))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB32))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB33))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB34)))
                trapz(csi,trapz(-y,(C_dk(:,:,k).*BB41))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB42))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB43))) , trapz(csi,trapz(-y,(C_dk(:,:,k).*BB44)))];
        end
        C = zeros(2*N,2*N,N-1);
        for k = 1:N-1
            C(:,:,k) = Extr_tr(:,:,k) * C_k(:,:,k) * Extr(:,:,k);
        end
        C = sum(C,3);
        CFF = C(1:2*N-2,1:2*N-2) + CFF_HLE;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xisp0(ex_indp) = xisp(ex_indp); xism0(ex_indp) = xism(ex_indp); xiT0(ex_indp) = xiT(ex_indp);
        xisp0(ex_indm) = xisp(ex_indm); xism0(ex_indm) = xism(ex_indm); xiT0(ex_indm) = xiT(ex_indm);
        status = 1;
    else
        status = 0;
    end
else
    switch(flag)
        case 'init'
            fprintf('start\n');
            tt = T(1);
            yy = X;
        case 'done'
            fprintf('done\n');
    end
    status = 0;
end
end
