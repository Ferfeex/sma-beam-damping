function status = myfun(T,X,flag,prm)
% OutputFcn of the ODE solver: this is where the SMA physics lives.
% After every accepted time step it
%   1. computes the strain field from the displacements (epsilon = B*X_k),
%   2. finds which integration points are inside a martensitic transformation
%      zone and solves the stress/martensite fractions there with Newton-Raphson,
%   3. rebuilds the tangent stiffness matrix KFF,
%   4. when a point leaves its transformation zone (stress reversal) it updates
%      the loss factor LF from the polynomial p and rebuilds the damping matrix CFF
%      (C = (eta/omega) K), then returns status=1 so that the solver restarts.
%
% Refactor note (2026): the original had 36 positional arguments; they are now
% fields of the struct prm and unpacked below with their ORIGINAL names.
% The numerical core (everything between the markers) is unchanged.
global KFF CFF tt yy xisp xism xiT xisp0 xism0 xiT0 indp_old indm_old sigma_old s p1 m1
global epsyylon xysp xysm xyT LF FF FFF
persistent last_print     % progress printing (not part of the physics)
Ea = prm.Ea; Emp = prm.Emp; Emm = prm.Emm; EmT = prm.EmT;
ssp = prm.ssp; sfp = prm.sfp; ssm = prm.ssm; sfm = prm.sfm;
eps_Lp = prm.eps_Lp; eps_Lm = prm.eps_Lm; b = prm.b;
Extr = prm.Extr; Extr_tr = prm.Extr_tr; csi = prm.csi; y = prm.y; B = prm.B; N = prm.N;
BB11 = prm.BB11; BB12 = prm.BB12; BB13 = prm.BB13; BB14 = prm.BB14;
BB21 = prm.BB21; BB22 = prm.BB22; BB23 = prm.BB23; BB24 = prm.BB24;
BB31 = prm.BB31; BB32 = prm.BB32; BB33 = prm.BB33; BB34 = prm.BB34;
BB41 = prm.BB41; BB42 = prm.BB42; BB43 = prm.BB43; BB44 = prm.BB44;
p = prm.p;
f_ref = prm.f_ref;          % reference frequency [Hz] in C = LF/(2*pi*f_ref)*K  (was hard-coded 10 / 14.5 / 18 ...)
KFF_ADD = prm.KFF_ADD;      % stiffness of the non-SMA layers (zeros for the pure-SMA beam)
CFF_ADD = prm.CFF_ADD;      % damping   of the non-SMA layers (zeros for the pure-SMA beam)
if nargin < 3 || isempty(flag)
    % VARIABLES UPDATE AND EVALUATION
    % External load (optional): trapezoidal pulse, ramp up in t_rise, hold until t_end, then zero.
    % (this block comes from Loss/Transient/myfun.m; with load.type='none' FFF stays zero)
    if prm.has_pulse
        if T <= prm.t_rise
            FFF = T*FF/prm.t_rise;
        elseif T > prm.t_rise && T < prm.t_end
            FFF = FF;
        else
            FFF = 0*FF;
        end
    end
    if prm.perf.hist, tt = [tt T]; yy = [yy X]; end   % history of the original (never read back; O(n^2) copy cost)
    % DISPLACEMENT FIELD
    XA = [X(2*N-1:end);zeros(2,1)];   % displacement part of the state (was X(31:60) for N=16)
    if prm.perf.vec
        epsilon = prm.Bfull * X(2*N-1:end);   % B*Extr precombined (identical result up to round-off)
    else
    X_k = zeros(4*(N-1),1);
    for k = 1:N-1
        X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * XA;
    end
    epsilon = B * X_k;
    end
    % FIRST SIGMA GUESS
    sigma_i = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon-eps_Lp*xisp-eps_Lm*xism);
    % NEWTON RAPHSON + FSOLVE (WHEN NEWTON DOES NOT CONVERGE)
    i = 0; % Iterations counter
    %%%%% LOCAL STATE ACCORDING TO MARTENSITIC PHASE DIAGRAM %%%%%%%%%%%%%%%%%%
    indp = find(sigma_i > ssp & sigma_i < sfp & sigma_i >= sigma_old);
    indm = find(sigma_i < ssm & sigma_i > sfm & sigma_i <= sigma_old);
    %%%%%%%%%%%%%%%%%%% UPDATING XI_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if prm.perf.sets
        % logical-mask version of the set algebra below (same index sets, no setdiff/union/unique overhead)
        n_ = numel(sigma_i);
        Mp = msk_(indp,n_); Mm = msk_(indm,n_);
        exp_ = msk_(indp_old,n_) & ~Mp;  exm_ = msk_(indm_old,n_) & ~Mm;
        ex_indp = find(exp_); ex_indm = find(exm_);
        P1 = msk_(p1,n_); M1 = msk_(m1,n_);
        M1(sigma_i > 0) = false;
        P1(sigma_i < 0) = false;
        P1 = P1 | exp_;  M1 = M1 | exm_;
        indp = find(Mp & ~exp_ & ~P1);  indm = find(Mm & ~exm_ & ~M1);
        p1 = find(P1); m1 = find(M1);
    else
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
    end
    
    indp_old = indp; indm_old = indm;
    
    sigmas = sigma_i;
    dF = zeros(size(xisp));   % was zeros(750,1): (N-1)*5*10 integration points
    Dxisp = zeros(size(xisp)); Dxism = zeros(size(xism)); DxiT = zeros(size(xiT));
    % NOTE (see KNOWN_ISSUES.md #1): the 2020 original had
    %     while ~isempty(indp) && ~isempty(indm)
    % i.e. Newton ran only while BOTH a tension and a compression zone were active and stopped as soon
    % as one of the two sets had converged. Corrected (2026, approved by the author) to ||: iterate while
    % ANY point still has to converge. The cap prm.newton_maxit protects against a non-converging point.
    while (~isempty(indp) || ~isempty(indm)) && i < prm.newton_maxit
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
        if prm.perf.sets
            indp = indp(abs(F(indp))>10^-6);  indm = indm(abs(F(indm))>10^-6);
        else
        indp = indp(ismember(indp,index));
        indm = indm(ismember(indm,index));
        end
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
    if prm.perf.vec
        KV = assemble_element_matrix(C_tot, prm.pf);
        KFF = KV + KFF_ADD;
    else
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
    KFF = K(1:2*N-2,1:2*N-2) + KFF_ADD;
    end
    if prm.perf.hist
    epsyylon = [epsyylon, epsilon];
    xysp = [xysp , xisp];
    xysm = [xysm , xism];
    xyT = [xyT, xiT];
    s = [s , sigma_i];
    end
    sigma_old = sigma_i;
    if prm.verbose
        fprintf('t = %f [s] ', T)
        fprintf('i = %i [s] \n', i)
    else
        if isempty(last_print) || T < last_print, last_print = -1; end
        if T - last_print >= 0.05
            last_print = T;
            fprintf('  t = %.3f s\n', T);
        end
    end
    if any(ex_indp)~=0 || any(ex_indm)~=0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW DAMPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LF(ex_indp) = polyval(p,epsilon(ex_indp));
        LF(ex_indm) = polyval(p,abs(epsilon(ex_indm)));
        C_d = C_tot .* LF / (2*pi*f_ref);
        if prm.perf.vec
            CV = assemble_element_matrix(C_d, prm.pf);
            CFF = CV + CFF_ADD;
        else
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
        CFF = C(1:2*N-2,1:2*N-2) + CFF_ADD;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xisp0(ex_indp) = xisp(ex_indp); xism0(ex_indp) = xism(ex_indp); xiT0(ex_indp) = xiT(ex_indp);
        xisp0(ex_indm) = xisp(ex_indm); xism0(ex_indm) = xism(ex_indm); xiT0(ex_indm) = xiT(ex_indm);
        status = 1;
    else
        status = 0;
        % restart the solver when the external load changes slope (from Loss/Transient/myfun.m)
        if prm.has_pulse && ((T >= prm.t_rise && T < prm.t_rise + prm.dt_restart) || (T >= prm.t_end && T < prm.t_end + prm.dt_restart))
            status = 1;
        end
    end
else
    switch(flag)
        case 'init'
            if prm.verbose, fprintf('start\n'); end
            tt = T(1);
            yy = X;
        case 'done'
            if prm.verbose, fprintf('done\n'); end
    end
    status = 0;
end
end

function z = msk_(ix,n)
z = false(n,1);  z(ix) = true;
end
