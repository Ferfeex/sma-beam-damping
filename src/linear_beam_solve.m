function res = linear_beam_solve(cfg)
% LINEAR_BEAM_SOLVE  linear (non-SMA) Euler-Bernoulli beam with Rayleigh damping.
%   This is the reference the SMA results are compared against ("HLE" = high-loss-elastic
%   linear material of same stiffness), and the linear transient cases (thesis 3.2/3.3.2).
%   Formerly: root FEM1D.m and Transient_linear/FEM1D.m (+ funforce.m).
%
%   cfg.mat.E, cfg.mat.rho, cfg.geom.*, cfg.damping.zeta (Rayleigh ratio on modes 1 and 2),
%   cfg.init (mode1_velocity | zero), cfg.load (none | pulse | ramp), cfg.time.Tf.
%
%   NOTE on the force: the stored original of Transient_linear/fun_ode_1D.m had no force term
%   (4 arguments, FEM1D.m called it with 6), so the form  M a = -K u - C v + F  with
%   F = FFF (ramped by funforce) is INFERRED, not recovered from the original.
global FF FFF
N = cfg.geom.N;
bm = build_beam(cfg, cfg.mat.E);
MFF = bm.MFF;  KFF = bm.KFF;  Extr_tr = bm.Extr_tr;

% Rayleigh damping: zeta on the first two modes
[eigenvectors, eigenvalues]=eig(MFF\KFF);
freq=sqrt(diag(eigenvalues))/2/pi;
[dummy,ind]=sort(freq);
h1=cfg.damping.zeta; h2=cfg.damping.zeta;
h=[h1;h2];
om1=2*pi*dummy(1);
om2=2*pi*dummy(2);
H=[1/(2*om1) om1/2;1/(2*om2) om2/2];
coef=H\h;
alpha=coef(1);
beta=coef(2);
CFF=alpha*MFF+beta*KFF;

% load
FF = zeros(2*N-2,1); FFF = 0*FF;
has_force = ~strcmp(cfg.load.type,'none');
if has_force
    F5 = zeros(4,1);
    F5(1) = cfg.load.F_element;
    F = Extr_tr(:,:,N/2) * F5 + Extr_tr(:,:,N/2+1) * F5;
    FF = F(1:2*N-2);
    FFF = 0*FF;
end
ld = cfg.load;

T0 = 0; Tf = cfg.time.Tf;
if strcmp(cfg.init.type,'mode1_velocity')
    X0 = initial_velocity(N, cfg.init.scale);
else
    X0 = zeros(4*N-4,1);
end

solver = cfg.solver.name;
if exist('OCTAVE_VERSION','builtin') ~= 0 && isfield(cfg.solver,'name_octave'), solver = cfg.solver.name_octave; end
solve = str2func(solver);
if has_force
    options = odeset('OutputFcn',@(T,X,flag) funforce(T,X,flag,ld),'InitialStep',cfg.solver.InitialStep);
else
    options = odeset();
end
rhs = @(T,X) rhs_linear(T,X,MFF,CFF,KFF,N);
tic_start = tic;
[T,X] = solve(rhs,[T0 Tf],X0,options);
TT = T;
XX = X;
T0 = T(end);
X0 = X(end,:)';
while T(end)<Tf-0.0001
    [T,X] = solve(rhs,[T0 Tf],X0,options);
    TT = [TT ; T(2:end)];
    XX = [XX ; X(2:end,:)];
    T0 = TT(end);
    X0 = X(end,:)';
end
col = 2*N-2+2*N/2-1;
if cfg.post.shift_to_first_max
    indice = find(XX(:,col) == max(XX(:,col)));
    indice = indice(1);
    TT = TT(indice:end)-TT(indice);
    XX = XX(indice:end,:);
end
res.TT = TT; res.XX = XX; res.col = col; res.midspan = XX(:,col);
res.runtime_s = toc(tic_start); res.alpha = alpha; res.beta = beta;
end

function xdot = rhs_linear(T,X,M,C,K,N)
% inferred force term (see header)
global FFF
xdot = [M\(-C*X(1:2*N-2) - K*X(2*N-1:end) + FFF)
    X(1:2*N-2)];
end

function status = funforce(T,~,flag,ld)
% Load history (formerly funforce.m): ramp from 0 to FF in t_rise, hold until t_end, then 0.
%   incremental pulse : t_rise = 0.033, t_end = 0.033  (load.type = 'ramp')
%   rectangular pulse : t_rise = 1e-4,  t_end = 0.022  (load.type = 'pulse')
% The solver is restarted (status=1) just after each change of slope.
global FF FFF
status = 0;
if ~isempty(flag), return; end
if T <= ld.t_rise
    FFF = T*FF/ld.t_rise;
elseif T > ld.t_rise && T < ld.t_end
    FFF = FF;
else
    FFF = 0*FF;
end
if (T > ld.t_rise && T < ld.t_rise + ld.dt_restart) || (T >= ld.t_end && T < ld.t_end + ld.dt_restart)
    status = 1;
end
end
