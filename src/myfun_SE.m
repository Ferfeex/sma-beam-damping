function status = myfun_SE(T,X,flag,prm)
% OutputFcn for the SUPER-ELASTIC, pre-stressed case (one martensite fraction xi_s^+ only).
% Same structure as myfun.m (strain -> Newton on stress/xi -> tangent K -> loss factor/damping
% update at stress reversal), with the austenite<->detwinned-martensite phase diagram of
% the super-elastic regime (T0 above Af). Taken from Loss/SE_pre-strain/myfun.m; the
% 42 positional arguments are now fields of prm (unpacked with their ORIGINAL names).
global KFF CFF tt yy xisp xisp0 indp_old indm_old sigma_old s p1 m1
global epsyylon xysp LF sigma0 epsilon0
persistent last_print
Ea = prm.Ea; Emp = prm.Emp; ssp = prm.ssp; sfp = prm.sfp; eps_Lp = prm.eps_Lp; b = prm.b;
Extr = prm.Extr; Extr_tr = prm.Extr_tr; csi = prm.csi; y = prm.y; B = prm.B; N = prm.N;
BB11 = prm.BB11; BB12 = prm.BB12; BB13 = prm.BB13; BB14 = prm.BB14;
BB21 = prm.BB21; BB22 = prm.BB22; BB23 = prm.BB23; BB24 = prm.BB24;
BB31 = prm.BB31; BB32 = prm.BB32; BB33 = prm.BB33; BB34 = prm.BB34;
BB41 = prm.BB41; BB42 = prm.BB42; BB43 = prm.BB43; BB44 = prm.BB44;
Temp0 = prm.Temp0; ssmart = prm.ssmart; sfmart = prm.sfmart; ssaus = prm.ssaus; sfaus = prm.sfaus;
As = prm.As; Af = prm.Af; CAp = prm.CAp; p = prm.p; eps0 = prm.eps0;
f_ref = prm.f_ref;                 % was hard-coded 21 [Hz]
LF_offset = prm.LF_offset;         % was hard-coded +0.015
sigma_rev = prm.sigma_rev;         % was hard-coded 120 [MPa]
if nargin < 3 || isempty(flag)
    % VARIABLES UPDATE AND EVALUATION
    if prm.perf.hist == 1, tt = [tt T]; yy = [yy X]; end   % history of the original (never read back; O(n^2) copy cost)
    % DISPLACEMENT FIELD
    XA = [X(2*N-1:end);zeros(2,1)];   % was X(31:60) for N=16
    if prm.perf.vec
        epsilon = prm.Bfull * X(2*N-1:end) + eps0;   % B*Extr precombined (identical result up to round-off)
    else
    X_k = zeros(4*(N-1),1);
    for k = 1:N-1
        X_k(4*(k-1)+1:4*k,1) = Extr(:,:,k) * XA;
    end
    epsilon = B * X_k + eps0;
    end
    % FIRST SIGMA GUESS
    sigma_i = sigma0 + (Ea + xisp*(Emp-Ea)) .* (epsilon-eps_Lp*xisp) - (Ea + xisp0*(Emp-Ea)) .* (epsilon0-eps_Lp*xisp0);
    % NEWTON RAPHSON + FSOLVE (WHEN NEWTON DOES NOT CONVERGE)
    i = 0; % Iterations counter
    %%%%% LOCAL STATE ACCORDING TO MARTENSITIC PHASE DIAGRAM %%%%%%%%%%%%%%%%%%
    indp = find(sigma_i > ssmart & sigma_i < sfmart & sigma_i >= sigma_old);
    indm = find(sigma_i < ssaus & sigma_i > sfaus & sigma_i <= sigma_old);
    %%%%%%%%%%%%%%%%%%% UPDATING XI_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if prm.perf.sets
        % logical-mask version of the set algebra below (same index sets, no setdiff/union/unique overhead)
        n_ = numel(sigma_i);
        Mp = msk_(indp,n_); Mm = msk_(indm,n_);
        exp_ = msk_(indp_old,n_) & ~Mp;  exm_ = msk_(indm_old,n_) & ~Mm;
        ex_indp = find(exp_); ex_indm = find(exm_);
        P1 = msk_(p1,n_); M1 = msk_(m1,n_);
        M1(sigma_i > ssaus) = false;
        P1(sigma_i < ssmart) = false;
        P1 = P1 | exp_;  M1 = M1 | exm_;
        indp = find(Mp & ~exp_ & ~P1);  indm = find(Mm & ~exm_ & ~M1);
        p1 = find(P1); m1 = find(M1);
    else
    ex_indp = setdiff(indp_old, indp); ex_indm = setdiff(indm_old, indm);
    am = find(sigma_i > ssaus);
    ap = find(sigma_i < ssmart);
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
    dF = zeros(size(xisp));
    Dxisp = zeros(size(xisp));
    % (same correction as in myfun.m: || instead of the original &&, see KNOWN_ISSUES.md #1)
    while (~isempty(indp) || ~isempty(indm)) && i < prm.newton_maxit
        % MARTENSITIC FRACTIONS
        xisp(indp) = (1-cos(pi*(sigma_i(indp)-ssmart)/(sfp-ssp)))/2 + xisp0(indp).*(1+cos(pi*(sigma_i(indp)-ssmart)/(sfp-ssp)))/2;
        xisp(indm) = xisp0(indm)/2 .* (cos(pi*(Temp0-As-sigma_i(indm)/CAp)/(Af-As))+1);
        % SIGMA FOR COMPARISON
        sigmas = sigma0 + (Ea + xisp*(Emp-Ea)) .* (epsilon-eps_Lp*xisp) - (Ea + xisp0*(Emp-Ea)) .* (epsilon0-eps_Lp*xisp0);
        % MARTENSITIC FRACTIONS DERIVATIVES
        Dxisp(indp) = (xisp0(indp)-1)/2 * pi/(ssp - sfp) .* sin(pi*(sigma_i(indp)-ssmart)/(sfp-ssp));
        Dxisp(indm) = xisp0(indm)/2 * pi/(CAp * (Af-As)) .* sin(pi*(Temp0 - As - sigma_i(indm)/CAp)/(Af - As));
        % YOUNG'S MODULUS DERIVATIVE
        dE = Dxisp*(Emp-Ea);
        % NEWTON'S RESIDUAL, EVALUATION & DERIVATIVE
        F = sigma_i-sigmas; % At those indeces ~=indp,indm this difference =0
        index = find(abs(F)>10^-12);
        dF(index) = 1 - (dE(index).*(epsilon(index)-eps_Lp*xisp(index)) + (Ea + xisp(index)*(Emp-Ea)).*(-eps_Lp*Dxisp(index)));
        % SIGMA UPDATE
        sigma_i(index) = sigma_i(index) - F(index)./dF(index);
        % UNION OF PHASE DIAGRAM & NON-CONVERGENCE INDECES
        if prm.perf.sets
            indp = indp(abs(F(indp))>10^-12);  indm = indm(abs(F(indm))>10^-12);
        else
        indp = indp(ismember(indp,index));
        indm = indm(ismember(indm,index));
        end
        i = i+1; % Counter update
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST NEWTON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigma_i = sigmas;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NEW YOUNG'S MODULUS
    E = Ea + xisp*(Emp-Ea);
    dE = Dxisp*(Emp-Ea);
    % TANGENT MATRIX & STIFFNESS MATRIX FOR K-TH ELEMENT
    C_tot = E ./ ( 1-(epsilon-eps_Lp*xisp).*dE + E.*(eps_Lp*Dxisp) );
    Ck = reshape(C_tot,length(y),length(csi),N-1);
    if prm.perf.vec
        KV = assemble_element_matrix(C_tot, prm.pf);
        KFF = KV;
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
    KFF = K(1:2*N-2,1:2*N-2);
    end
    if prm.perf.hist == 1          % full history, one column per step (O(n^2) copying: only if cfg.post.keep_field_history)
    epsyylon = [epsyylon, epsilon];
    xysp = [xysp , xisp];
    s = [s , sigma_i];
    elseif prm.perf.hist == 2      % only the latest column (all that the field plots use)
        epsyylon = epsilon;
        xysp = xisp;
        s = sigma_i;
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
        LF(ex_indp) = polyval(p,sigma_i(ex_indp))+LF_offset;
        LF(ex_indm) = polyval(p,sigma_rev - sigma_i(ex_indm))+LF_offset;
        C_d = C_tot .* LF / (2*pi*f_ref);
        if prm.perf.vec
            CV = assemble_element_matrix(C_d, prm.pf);
            CFF = CV;
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
        CFF = C(1:2*N-2,1:2*N-2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xisp0(ex_indp) = xisp(ex_indp);
        xisp0(ex_indm) = xisp(ex_indm);
        sigma0(ex_indp) = sigma_i(ex_indp); epsilon0(ex_indp) = epsilon(ex_indp);
        sigma0(ex_indm) = sigma_i(ex_indm); epsilon0(ex_indm) = epsilon(ex_indm);
        status = 1;
    else
        status = 0;
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
