function res = run_case(c, varargin)
% RUN_CASE  single entry point: run one simulation case and (optionally) plot / save it.
%
%   run_case('base')                      free vibration of the SMA beam (thesis base case)
%   run_case('composite')                 GFRP/SMA/GFRP sandwich beam
%   run_case('hle_linear')                linear reference beam (E = 26 GPa, Rayleigh damping)
%   run_case('pulse')                     rectangular force pulse (SLOW, see README)
%   run_case('prestrain_se')              super-elastic pre-stressed beam (SLOW)
%   run_case('linear_pulse_incremental'), run_case('linear_pulse_rectangular')
%   run_case(cfg)                         cfg = struct (e.g. filled by a future menu)
%
%   Name/value options:
%     'Tf'          override the duration [s]        (e.g. 'Tf',0.1 for a quick smoke test)
%     'compare_hle' true/false: also run the linear reference and overlay it (default false)
%     'figures'     'save' (default, PNG in outdir) | 'show' | 'none'
%     'outdir'      folder for figures/results (default <repo>/results, git-ignored)
%     'save_mat'    true/false: save TT, XX, cfg to outdir (default true)
%     'newton_or'   true = use || in the Newton loop condition (KNOWN_ISSUES.md #1)
%     'verbose'     true = original per-step printout
%
%   Returns res with fields TT, XX, midspan, runtime_s, cfg, peaks (time, height).
setup_paths();
opt = struct('Tf',[],'compare_hle',false,'figures','save','outdir','', ...
    'save_mat',true,'newton_or',[],'verbose',[]);
for k = 1:2:numel(varargin)
    if ~isfield(opt,varargin{k}), error('run_case: unknown option "%s"', varargin{k}); end
    opt.(varargin{k}) = varargin{k+1};
end
root = fileparts(mfilename('fullpath'));
if isempty(opt.outdir), opt.outdir = fullfile(root,'results'); end

if ischar(c)
    cfg = feval(['case_' c]);
else
    cfg = c;
end
if ~isempty(opt.Tf),        cfg.time.Tf = opt.Tf; end
if ~isempty(opt.newton_or), cfg.numerics.newton_or = opt.newton_or; end
if ~isempty(opt.verbose),   cfg.numerics.verbose = opt.verbose; end

if ~strcmp(cfg.model,'linear')
    d = load(fullfile(root,'data',cfg.p_file));   % polynomial p of the loss factor
    cfg.p = d.p;
end

fprintf('== run_case: %s ==\n%s\n', cfg.name, cfg.description);
if strcmp(cfg.model,'linear')
    res = linear_beam_solve(cfg);
else
    res = sma_beam_solve(cfg);
end
[tp,xp] = sma_peaks(res.TT,res.midspan);
res.peaks = [tp xp];
res.cfg = cfg;
fprintf('runtime %.1f s, %d time steps, %d solver restarts\n', res.runtime_s, numel(res.TT), getfield_or(res,'nrestart',0));
fprintf('mid-span peaks [mm] (t[s] : height): ');
fprintf('%.3f:%.2f  ', res.peaks(1:min(8,end),:)');
fprintf('\n');

ref = [];
if opt.compare_hle
    cr = case_hle_linear(); cr.time.Tf = cfg.time.Tf;
    ref = linear_beam_solve(cr);
    res.reference = ref;
end

if ~strcmp(opt.figures,'none')
    if strcmp(opt.figures,'save'), f = figure('visible','off'); else f = figure; end
    plot(res.TT, res.midspan, 'LineWidth', 1.5); hold on
    if ~isempty(ref), plot(ref.TT, ref.midspan, 'LineWidth', 1.5); end
    grid on; xlabel('Time t [s]'); ylabel('Displacement \delta [mm]');
    title(sprintf('Mid-span vertical displacement - %s', strrep(cfg.name,'_','\_')));
    if ~isempty(ref), legend('SMA (loss-factor damping)','linear reference, E = 26 GPa, Rayleigh 1.3 %'); end
    if strcmp(opt.figures,'save')
        if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
        fn = fullfile(opt.outdir, [cfg.name '_midspan.png']);
        print(f,'-dpng','-r110',fn); close(f);
        fprintf('figure: %s\n', fn);
    end
end
if opt.save_mat
    if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
    TT = res.TT; XX = res.XX; midspan = res.midspan; peaks = res.peaks; %#ok<NASGU>
    save('-v7', fullfile(opt.outdir,[cfg.name '_result.mat']), 'TT','XX','midspan','peaks','cfg');
end
end

function v = getfield_or(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
