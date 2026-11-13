function status = myfun_SE(T,X,flag,prm)
% OutputFcn for the SUPER-ELASTIC, pre-stressed case (one martensite fraction xi_s^+ only).
% Same structure as myfun.m (strain -> Newton on stress/xi -> tangent K -> loss factor/damping
% update at stress reversal), with the austenite<->detwinned-martensite phase diagram of
% the super-elastic regime (T0 above Af). The parameters come in the struct prm (see sma_beam_solve.m).
global KFF CFF xisp xisp0 indp_old indm_old sigma_old s p1 m1
global epsyylon xysp LF sigma0 epsilon0
global snap_n snap_ns snap_gi snap_grid snap_T snap_P snap_M snap_Tm
persistent last_print
Ea = prm.Ea; Emp = prm.Emp; ssp = prm.ssp; sfp = prm.sfp; eps_Lp = prm.eps_Lp;
N = prm.N;
Temp0 = prm.Temp0; ssmart = prm.ssmart; sfmart = prm.sfmart; ssaus = prm.ssaus; sfaus = prm.sfaus;
As = prm.As; Af = prm.Af; CAp = prm.CAp; p = prm.p; eps0 = prm.eps0;
f_ref = prm.f_ref;                 % reference frequency [Hz]
LF_offset = prm.LF_offset;         % offset added to the loss factor
sigma_rev = prm.sigma_rev;         % stress [MPa] used to mirror the unloading branch of the loss-factor curve
if nargin < 3 || isempty(flag)
    % VARIABLES UPDATE AND EVALUATION
    % DISPLACEMENT FIELD
    epsilon = prm.Bfull * X(2*N-1:end) + eps0;   % strain at the integration points (B*Extr precombined) plus the pre-strain
    % FIRST SIGMA GUESS
    sigma_i = sigma0 + (Ea + xisp*(Emp-Ea)) .* (epsilon-eps_Lp*xisp) - (Ea + xisp0*(Emp-Ea)) .* (epsilon0-eps_Lp*xisp0);
    % NEWTON RAPHSON + FSOLVE (WHEN NEWTON DOES NOT CONVERGE)
    i = 0; % Iterations counter
    %%%%% LOCAL STATE ACCORDING TO MARTENSITIC PHASE DIAGRAM %%%%%%%%%%%%%%%%%%
    indp = find(sigma_i > ssmart & sigma_i < sfmart & sigma_i >= sigma_old);
    indm = find(sigma_i < ssaus & sigma_i > sfaus & sigma_i <= sigma_old);
    %%%%%%%%%%%%%%%%%%% UPDATING XI_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % index sets of the transformation zones (logical masks)
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

    indp_old = indp; indm_old = indm;

    sigmas = sigma_i;
    dF = zeros(size(xisp));
    Dxisp = zeros(size(xisp));
    % Newton iterations run while ANY point still has to converge (as in myfun.m)
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
            indp = indp(abs(F(indp))>10^-12);  indm = indm(abs(F(indm))>10^-12);
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
    KV = assemble_element_matrix(C_tot, prm.pf);
    KFF = KV;
    if prm.histmode == 1          % full history, one column per step (O(n^2) copying: only if cfg.post.keep_field_history)
        epsyylon = [epsyylon, epsilon];
        xysp = [xysp , xisp];
        s = [s , sigma_i];
    elseif prm.histmode == 2      % only the latest column (all that the field plots use)
        epsyylon = epsilon;
        xysp = xisp;
        s = sigma_i;
    end
    if snap_gi <= snap_n && T >= snap_grid(snap_gi)      % optional snapshot of the martensite fractions (cfg.post.n_snapshots)
        snap_ns = snap_ns + 1;
        snap_T(snap_ns) = T;  snap_P(:,snap_ns) = xisp;
        while snap_gi <= snap_n && snap_grid(snap_gi) <= T, snap_gi = snap_gi + 1; end
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
        CV = assemble_element_matrix(C_d, prm.pf);
        CFF = CV;
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
        case 'done'
            if prm.verbose, fprintf('done\n'); end
    end
    status = 0;
end
end

function z = msk_(ix,n)
z = false(n,1);  z(ix) = true;
end
