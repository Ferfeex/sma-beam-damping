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
% The parameters come in the struct prm (built in sma_beam_solve.m). The state of the martensite fractions
% and of the damping is kept in globals, which sma_beam_solve resets at the start of every solve.
global KFF CFF xisp xism xiT xisp0 xism0 xiT0 indp_old indm_old sigma_old s p1 m1
global epsyylon xysp xysm xyT LF FF FFF
global snap_n snap_ns snap_gi snap_grid snap_T snap_P snap_M snap_Tm
persistent last_print     % progress printing (not part of the physics)
Ea = prm.Ea; Emp = prm.Emp; Emm = prm.Emm; EmT = prm.EmT;
ssp = prm.ssp; sfp = prm.sfp; ssm = prm.ssm; sfm = prm.sfm;
eps_Lp = prm.eps_Lp; eps_Lm = prm.eps_Lm;
N = prm.N;
p = prm.p;
f_ref = prm.f_ref;          % reference frequency [Hz] in C = LF/(2*pi*f_ref)*K
KFF_ADD = prm.KFF_ADD;      % stiffness of the non-SMA layers (zeros for the pure-SMA beam)
CFF_ADD = prm.CFF_ADD;      % damping   of the non-SMA layers (zeros for the pure-SMA beam)
if nargin < 3 || isempty(flag)
    % VARIABLES UPDATE AND EVALUATION
    % External load (optional): trapezoidal pulse, ramp up in t_rise, hold until t_end, then zero.
    % (with load.type = 'none' FFF stays zero)
    if prm.has_pulse
        if T <= prm.t_rise
            FFF = T*FF/prm.t_rise;
        elseif T > prm.t_rise && T < prm.t_end
            FFF = FF;
        else
            FFF = 0*FF;
        end
    end
    % DISPLACEMENT FIELD
    epsilon = prm.Bfull * X(2*N-1:end);   % strain at the integration points from the displacement part of the state (B*Extr precombined)
    % FIRST SIGMA GUESS
    sigma_i = (Ea + xisp*(Emp-Ea) + xism*(Emm-Ea) + xiT*(EmT-Ea)) .* (epsilon-eps_Lp*xisp-eps_Lm*xism);
    % NEWTON RAPHSON + FSOLVE (WHEN NEWTON DOES NOT CONVERGE)
    i = 0; % Iterations counter
    %%%%% LOCAL STATE ACCORDING TO MARTENSITIC PHASE DIAGRAM %%%%%%%%%%%%%%%%%%
    indp = find(sigma_i > ssp & sigma_i < sfp & sigma_i >= sigma_old);
    indm = find(sigma_i < ssm & sigma_i > sfm & sigma_i <= sigma_old);
    %%%%%%%%%%%%%%%%%%% UPDATING XI_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % index sets of the transformation zones (logical masks)
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

    indp_old = indp; indm_old = indm;

    sigmas = sigma_i;
    dF = zeros(size(xisp));
    Dxisp = zeros(size(xisp)); Dxism = zeros(size(xism)); DxiT = zeros(size(xiT));
    % Newton iterations run while ANY point still has to converge; prm.newton_maxit caps them.
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
            indp = indp(abs(F(indp))>10^-6);  indm = indm(abs(F(indm))>10^-6);
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
    KV = assemble_element_matrix(C_tot, prm.pf);
    KFF = KV + KFF_ADD;
    if prm.histmode == 1          % full history, one column per step (O(n^2) copying: only if cfg.post.keep_field_history)
        epsyylon = [epsyylon, epsilon];
        xysp = [xysp , xisp];
        xysm = [xysm , xism];
        xyT = [xyT, xiT];
        s = [s , sigma_i];
    elseif prm.histmode == 2      % only the latest column (all that the field plots use)
        epsyylon = epsilon;
        xysp = xisp;
        xysm = xism;
        xyT = xiT;
        s = sigma_i;
    end
    if snap_gi <= snap_n && T >= snap_grid(snap_gi)      % optional snapshot of the martensite fractions (cfg.post.n_snapshots)
        snap_ns = snap_ns + 1;
        snap_T(snap_ns) = T;  snap_P(:,snap_ns) = xisp;
        snap_M(:,snap_ns) = xism;  snap_Tm(:,snap_ns) = xiT;
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
        LF(ex_indp) = polyval(p,epsilon(ex_indp));
        LF(ex_indm) = polyval(p,abs(epsilon(ex_indm)));
        C_d = C_tot .* LF / (2*pi*f_ref);
        CV = assemble_element_matrix(C_d, prm.pf);
        CFF = CV + CFF_ADD;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xisp0(ex_indp) = xisp(ex_indp); xism0(ex_indp) = xism(ex_indp); xiT0(ex_indp) = xiT(ex_indp);
        xisp0(ex_indm) = xisp(ex_indm); xism0(ex_indm) = xism(ex_indm); xiT0(ex_indm) = xiT(ex_indm);
        status = 1;
    else
        status = 0;
        % restart the solver when the external load changes slope
        if prm.has_pulse && ((T >= prm.t_rise && T < prm.t_rise + prm.dt_restart) || (T >= prm.t_end && T < prm.t_end + prm.dt_restart))
            status = 1;
        end
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
