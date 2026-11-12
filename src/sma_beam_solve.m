function res = sma_beam_solve(cfg)
% SMA_BEAM_SOLVE  Transient vibration of a simply supported Cu-Al-Zn-Mn SMA beam
%   (final method of the thesis, formerly Loss/FEM1D.m and its 5 copies).
%
%   The beam is solved as a LINEAR ODE system  M a + C v + K u = F  that is
%   advanced with ode23t (ode23s in Octave); the SMA nonlinearity enters through
%   the OutputFcn (myfun.m / myfun_SE.m), which after every step updates the
%   stiffness K(strain, martensite fractions) and, at every stress reversal, the
%   damping C = (eta/omega) K with the loss factor eta(strain) from the Brinson fit p.
%
%   res = sma_beam_solve(cfg)   cfg: struct from a cases/case_*.m file (see case_base.m)
%   res.TT, res.XX   time [s] and state [velocities(30) displacements(30)]
%   res.midspan      mid-span vertical displacement [mm]
%   res.fields       (optional) strain / martensite fractions / stress at every step
%
% Variable names (N, Ea, Emp, ..., KFF, CFF, LF, xisp, ...) are the original ones.

% ---- globals used by the OutputFcn (the original relied on "clear all" to reset them) ----
global KFF CFF LF FF FFF
global xisp xism xiT xisp0 xism0 xiT0 indp_old indm_old sigma_old s p1 m1 tt yy
global epsyylon xysp xysm xyT sigma0 epsilon0
p1 = []; m1 = []; tt = []; yy = [];
FF = []; FFF = 0;
sigma0 = []; epsilon0 = [];

model = cfg.model;
N = cfg.geom.N;
b = cfg.geom.b;
mat = cfg.mat;
p = cfg.p;                       % polynomial fit of the loss factor (Brinson/Tension_compression_parametric.m)

% ---- SMA properties ----
Ea = mat.Ea;  Emp = mat.Emp;  ssp = mat.ssp;  sfp = mat.sfp;  eps_Lp = mat.eps_Lp;
if strcmp(model,'two_phase')
    Emm = mat.Emm;  EmT = mat.EmT;  ssm = mat.ssm;  sfm = mat.sfm;  eps_Lm = mat.eps_Lm;
    E_K = mat.E_K;               % modulus used in the initial (linear) stiffness matrix
else
    E_K = Ea;
end

bm = build_beam(cfg, E_K);
Extr = bm.Extr; Extr_tr = bm.Extr_tr; csi = bm.csi; y = bm.y; B = bm.B;
MFF = bm.MFF;  KFF = bm.KFF;
nip = (N-1)*numel(csi)*numel(y);        % integration points (750 for N=16)

% ---- optional extra layers (composite: GFRP core + SMA layers) ----
MFF_ADD = zeros(size(MFF));  KFF_ADD = zeros(size(MFF));  CFF_ADD = zeros(size(MFF));
if isfield(cfg,'hle') && ~isempty(cfg.hle)
    if isfield(cfg.hle,'matrices_file') && ~isempty(cfg.hle.matrices_file)
        hle = load(cfg.hle.matrices_file);     % e.g. the 2020 Matrices.mat (used by the tests)
    else
        hle = composite_matrices(cfg.hle);
    end
    MFF_ADD = hle.MFF_HLE;  KFF_ADD = hle.KFF_HLE;  CFF_ADD = hle.CFF_HLE;
end
MFF = MFF + MFF_ADD;
KFF_SMA = KFF;
KFF = KFF_SMA + KFF_ADD;

% ---- external load (optional) ----
has_pulse = strcmp(cfg.load.type,'pulse');
if has_pulse
    % VERTICAL LOAD on the two elements around mid-span (F_element = nodal force per element)
    F5 = zeros(4,1);
    F5(1) = cfg.load.F_element;
    F = Extr_tr(:,:,N/2) * F5 + Extr_tr(:,:,N/2+1) * F5;
    FF = F(1:2*N-2);
    FFF = 0*FF;
end

% ---- initial damping from the loss factor at a reference strain (or stress) ----
f_ref = cfg.damping.f_ref;
LF_max = polyval(p,cfg.damping.x_ref);
CFF = KFF_SMA * LF_max / (2*pi*f_ref) + CFF_ADD;
LF = LF_max * ones(nip,1);

% ---- time interval and initial conditions ----
T0 = 0; Tf = cfg.time.Tf;
if strcmp(cfg.init.type,'mode1_velocity')
    X0 = initial_velocity(N, cfg.init.scale);
else
    X0 = zeros(4*N-4,1);
end

% ---- initial state of the martensite fractions ----
if strcmp(model,'two_phase')
    sigma_old = zeros(nip,1); xisp0 = zeros(nip,1); s = sigma_old;
    xism0 = zeros(nip,1); xiT0 = ones(nip,1); epsilon = zeros(nip,1);
    xisp = xisp0; xism = xism0; xiT = xiT0; indp_old = []; indm_old = [];
    xysp = xisp; xysm = xism; xyT = xiT; epsyylon = zeros(nip,1);
    eps0 = 0;
else
    % super-elastic case with pre-stress sigma0 (stress-free reference: austenite)
    Mf = mat.Mf; Ms = mat.Ms; As = mat.As; Af = mat.Af; CMp = mat.CMp; CAp = mat.CAp;
    Temp0 = cfg.init.Temp0;
    ssmart = CMp*(Temp0-Ms)+ssp;
    sfmart = CMp*(Temp0-Ms)+sfp;
    ssaus = CAp*(Temp0-As);
    sfaus = CAp*(Temp0-Af);
    E = Ea;
    sigma0 = cfg.init.sigma0 * ones(nip,1);
    epsilon0 = sigma0 ./ E;
    eps0 = cfg.init.sigma0/E;
    sigma_old = sigma0; xisp0 = zeros(nip,1); s = sigma_old; epsilon = epsilon0;
    xisp = xisp0; indp_old = []; indm_old = [];
    xysp = xisp; epsyylon = epsilon;
    xism = []; xiT = []; xysm = []; xyT = [];
end

% ---- PERFORMANCE OPTIONS (branch perf-experiments) ----
pfo = struct('vec',1,'hist',0,'fast_rhs',1,'jac',1,'sets',1);   % defaults = fast path (see PERF_NOTES.md); cfg.perf.<flag> = 0 restores the original code path
if isfield(cfg,'perf'), fn_ = fieldnames(cfg.perf); for k_ = 1:numel(fn_), pfo.(fn_{k_}) = cfg.perf.(fn_{k_}); end, end
pf = [];
if pfo.vec
    npt = numel(csi)*numel(y);
    wy = trapz(-y, eye(numel(y)));  wc = trapz(csi, eye(numel(csi)));
    W = wy(:)*wc;                                   % weights of trapz(csi,trapz(-y,.))
    BBc = {bm.BB11,bm.BB12,bm.BB13,bm.BB14, bm.BB21,bm.BB22,bm.BB23,bm.BB24, ...
           bm.BB31,bm.BB32,bm.BB33,bm.BB34, bm.BB41,bm.BB42,bm.BB43,bm.BB44};
    WB = zeros(npt,16);                             % column a+4*(c-1) <-> BB_ac
    for a_ = 1:4, for c_ = 1:4, WB(:,a_+4*(c_-1)) = W(:).*BBc{4*(a_-1)+c_}(:); end, end
    I_ = zeros(16*(N-1),1); J_ = I_;
    for k_ = 1:N-1
        g = zeros(4,1);
        for a_ = 1:4, g(a_) = find(Extr(a_,:,k_)); end
        [aa,cc] = ndgrid(1:4,1:4);
        I_((k_-1)*16+(1:16)) = g(aa(:));  J_((k_-1)*16+(1:16)) = g(cc(:));
    end
    Eall = zeros(4*(N-1),2*N-2);
    for k_ = 1:N-1, Eall(4*(k_-1)+1:4*k_,:) = Extr(:,1:2*N-2,k_); end   % rows of Extr, free dofs only
    pf = struct('b',b,'npt',npt,'nel',N-1,'WB',WB,'I',I_,'J',J_,'ndof',2*N,'nfree',2*N-2);
    Bfull = B*Eall;
else
    Bfull = [];
end

% ---- parameters handed to the OutputFcn ----
prm = struct('N',N,'b',b,'Extr',Extr,'Extr_tr',Extr_tr,'csi',csi,'y',y,'B',B, ...
    'BB11',bm.BB11,'BB12',bm.BB12,'BB13',bm.BB13,'BB14',bm.BB14, ...
    'BB21',bm.BB21,'BB22',bm.BB22,'BB23',bm.BB23,'BB24',bm.BB24, ...
    'BB31',bm.BB31,'BB32',bm.BB32,'BB33',bm.BB33,'BB34',bm.BB34, ...
    'BB41',bm.BB41,'BB42',bm.BB42,'BB43',bm.BB43,'BB44',bm.BB44, ...
    'Ea',Ea,'Emp',Emp,'ssp',ssp,'sfp',sfp,'eps_Lp',eps_Lp,'p',p,'f_ref',f_ref, ...
    'newton_maxit',cfg.numerics.newton_maxit, ...
    'verbose',cfg.numerics.verbose, ...
    'perf',struct('vec',pfo.vec,'sets',pfo.sets,'hist',pfo.hist || cfg.post.keep_fields),'pf',pf,'Bfull',Bfull, ...
    'has_pulse',has_pulse,'t_rise',0,'t_end',0,'dt_restart',1e-5);
if has_pulse
    prm.t_rise = cfg.load.t_rise;  prm.t_end = cfg.load.t_end;
    if isfield(cfg.load,'dt_restart'), prm.dt_restart = cfg.load.dt_restart; end
end
if strcmp(model,'two_phase')
    prm.Emm = Emm; prm.EmT = EmT; prm.ssm = ssm; prm.sfm = sfm; prm.eps_Lm = eps_Lm;
    prm.KFF_ADD = KFF_ADD; prm.CFF_ADD = CFF_ADD;
    outfun = @(T,X,flag) myfun(T,X,flag,prm);
else
    prm.Temp0 = Temp0; prm.ssmart = ssmart; prm.sfmart = sfmart; prm.ssaus = ssaus; prm.sfaus = sfaus;
    prm.As = As; prm.Af = Af; prm.CAp = CAp; prm.eps0 = eps0;
    prm.LF_offset = cfg.damping.LF_offset; prm.sigma_rev = cfg.damping.sigma_rev;
    outfun = @(T,X,flag) myfun_SE(T,X,flag,prm);
end

% ---- ODE solver options ----
sv = cfg.solver;
solver = sv.name;
if is_octave() && isfield(sv,'name_octave'), solver = sv.name_octave; end
extra = {};
if ~isempty(sv.RelTol), extra = [extra, {'RelTol',sv.RelTol,'AbsTol',sv.AbsTol}]; end
options2 = odeset('OutputFcn',outfun,extra{:});
if ~isempty(sv.InitialStep)
    options = odeset('OutputFcn',outfun,'InitialStep',sv.InitialStep,extra{:});
else
    options = options2;
end

% ---- TIME INTEGRATION (restarts whenever the OutputFcn returns status=1) ----
if pfo.fast_rhs
    Minv = inv(MFF);
    rhs = @(T,X) fun_ode_fast(T,X,Minv,N);
else
    rhs = @(T,X) fun_ode_1D(T,X,MFF,N);
end
if pfo.jac
    Minv = inv(MFF);
    jopt = {'Jacobian', @(T,X) jac_ode_fast(T,X,Minv,N)};
    options2 = odeset(options2, jopt{:});
    options  = odeset(options,  jopt{:});
end
if isfield(pfo,'solver') && ~isempty(pfo.solver), solver = pfo.solver; end
if isfield(pfo,'abstol') && ~isempty(pfo.abstol)
    options2 = odeset(options2,'AbsTol',pfo.abstol); options = odeset(options,'AbsTol',pfo.abstol);
end
if isfield(pfo,'reltol') && ~isempty(pfo.reltol)
    options2 = odeset(options2,'RelTol',pfo.reltol); options = odeset(options,'RelTol',pfo.reltol);
end
solve = str2func(solver);
tic_start = tic;
[T,X] = solve(rhs,[T0 Tf],X0,options);
TT = T;
XX = X;
T0 = T(end);
X0 = X(end,:)';
nrestart = 0;
while T(end) < Tf-0.0001
    [T,X] = solve(rhs,[T0 Tf],X0,options2);
    TT = [TT ; T(2:end)];
    XX = [XX ; X(2:end,:)];
    T0 = TT(end);
    X0 = X(end,:)';
    nrestart = nrestart + 1;
    if nrestart > cfg.solver.max_restarts
        error('sma_beam_solve: too many solver restarts (%d), stopping.', nrestart);
    end
end

% ---- MID-SPAN VERTICAL DISPLACEMENT ----
col = 2*N-2+2*N/2-1;                      % column of the mid-span displacement in X
if cfg.post.shift_to_first_max
    % the time axis starts at the first (largest) mid-span maximum, as in the original
    indice = find(XX(:,col) == max(XX(:,col)));
    indice = indice(1);
    TT = TT(indice:end)-TT(indice);
    XX = XX(indice:end,:);
    if pfo.hist || cfg.post.keep_fields      % field history exists only if it was recorded
        epsyylon = epsyylon(:,indice:end);
        xysp = xysp(:,indice:end);
        if strcmp(model,'two_phase')
            xysm = xysm(:,indice:end);
            xyT = xyT(:,indice:end);
        end
        s = s(:,indice:end);
    end
end

res.TT = TT;
res.XX = XX;
res.col = col;
res.midspan = XX(:,col);
res.runtime_s = toc(tic_start);
res.nrestart = nrestart;
if cfg.post.keep_fields
    res.fields = struct('epsyylon',epsyylon,'xysp',xysp,'xysm',xysm,'xyT',xyT,'s',s);
end
end

function tf = is_octave()
tf = exist('OCTAVE_VERSION','builtin') ~= 0;
end
