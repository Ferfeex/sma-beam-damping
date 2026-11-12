function res = run_case(c, varargin)
% RUN_CASE  single entry point: run one simulation case and (optionally) plot / save it.
%
%   run_case                              NO arguments: opens the interactive start window (sma_gui.m):
%                                         choose the case, edit the loss-factor polynomial p, tick
%                                         "compare with linear beam", press Run
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
%     'figures'     'save' (default, PNG in outdir, no window) | 'show' (windows only) |
%                   'both' (windows AND PNG, what the start window uses by default) | 'none'
%     'p'           loss-factor polynomial to use for THIS run (vector, polyval order) instead of
%                   data/p.mat (or p_SE.mat); data files are never modified
%     'outdir'      folder for figures/results (default <repo>/results, git-ignored)
%     'save_mat'    true/false: save TT, XX, cfg to outdir (default true)
%     'verbose'     true = original per-step printout
%     'plot_fields' true = also plot, at the last time step, the martensite fractions (src/meshare.m, <case>_fields.png)
%                   and the deformed beam + strain map (src/strain_figure.m, <case>_strain.png)
%     'animate'     false (default) | true = play the beam vibration in a window (src/animate_beam.m) |
%                   'gif' = write <outdir>/<case>_animation.gif (no window) | 'both' = window AND GIF.
%                   Slow motion, 25 fps, deflection exaggerated; works for every case.
%     'animate_opts' cell of extra name/value options for animate_beam, e.g. {'speed',0.05,'strain',false}
%
%   Returns res with fields TT, XX, midspan (deflection at x = L/2, Hermite-interpolated), midspan_node
%   (old value: node 9, x = res.x_midspan_node), runtime_s, cfg, peaks (time, height of midspan).
setup_paths();
if nargin == 0
    fig = sma_gui();                      % interactive start window (or text menu without display)
    if nargout > 0, res = fig; end
    return
end
opt = struct('Tf',[],'compare_hle',false,'figures','save','outdir','', ...
    'save_mat',true,'verbose',[],'plot_fields',false,'p',[],'animate',false,'animate_opts',{{}});
for k = 1:2:numel(varargin)
    if ~isfield(opt,varargin{k}), error('run_case: unknown option "%s"', varargin{k}); end
    opt.(varargin{k}) = varargin{k+1};
end
root = fileparts(mfilename('fullpath'));
if isempty(opt.outdir), opt.outdir = fullfile(root,'results'); end

if ~any(strcmp(opt.figures,{'save','show','both','none'}))
    error('run_case: option ''figures'' must be ''save'', ''show'', ''both'' or ''none''');
end
if ischar(c)
    cfg = feval(['case_' c]);
else
    cfg = c;
end
if ~isempty(opt.Tf),        cfg.time.Tf = opt.Tf; end
if ~isempty(opt.verbose),   cfg.numerics.verbose = opt.verbose; end
if opt.plot_fields,         cfg.post.keep_fields = true; end

if ~strcmp(cfg.model,'linear')
    pf = cfg.p_file;
    if ~exist(pf,'file'), pf = fullfile(root,'data',pf); end
    d = load(pf);                                  % polynomial p of the loss factor
    cfg.p = d.p;
    cfg.p_modified = false;
    if ~isempty(opt.p)
        if ~isnumeric(opt.p) || ~all(isfinite(opt.p(:))) || isempty(opt.p)
            error('run_case: option ''p'' must be a numeric vector of finite coefficients');
        end
        cfg.p = opt.p(:).';                        % this run only: nothing is written to data/
        cfg.p_modified = ~isequal(cfg.p, d.p(:).');
    end
else
    cfg.p_modified = false;
end
tag = '';  if cfg.p_modified, tag = '_pmod'; end   % do not overwrite the default-p files

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

show = any(strcmp(opt.figures,{'show','both'}));
dosave = any(strcmp(opt.figures,{'save','both'}));
if ~strcmp(opt.figures,'none')
    if show, f = figure; else, f = figure('visible','off'); end
    plot(res.TT, res.midspan, 'LineWidth', 1.5); hold on   % deflection at x = L/2
    if ~isempty(ref), plot(ref.TT, ref.midspan, 'LineWidth', 1.5); end
    grid on; xlabel('Time t [s]'); ylabel('Displacement \delta [mm]');
    ttl = sprintf('Mid-span vertical displacement (x = %g mm) - %s', cfg.geom.L/2, strrep(cfg.name,'_','\_'));
    if cfg.p_modified, ttl = [ttl ' (modified p)']; end
    title(ttl);
    if ~isempty(ref), legend(sprintf('SMA (loss-factor damping), x = %g mm',cfg.geom.L/2),sprintf('linear reference, E = 26 GPa, Rayleigh 1.3 %%, x = %g mm',cfg.geom.L/2)); end
    if dosave
        if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
        fn = fullfile(opt.outdir, [cfg.name tag '_midspan.png']);
        print(f,'-dpng','-r110',fn);
        if ~show, close(f); end
        fprintf('figure: %s\n', fn);
    end
    if show, drawnow; end
end
if opt.plot_fields && isfield(res,'fields') && strcmp(cfg.model,'two_phase') && ~strcmp(opt.figures,'none')
    ff = res.fields;
    if ~exist(opt.outdir,'dir') && dosave, mkdir(opt.outdir); end
    % 1) martensite fractions over the beam at the last time step (thesis 'Simulazione_lunga')
    if show, f = figure; else, f = figure('visible','off'); end
    meshare(ff.epsyylon,ff.xysp,ff.xysm,ff.xyT,res.TT,cfg.geom.thick,cfg.geom.L,f);
    if dosave
        fn = fullfile(opt.outdir, [cfg.name tag '_fields.png']);
        print(f,'-dpng','-r110',fn);
        fprintf('figure: %s\n', fn);
    end
    if ~show, close(f); else, drawnow; end
    % 2) thesis-style deformed beam + strain map at the last time step
    if show, f = figure; else, f = figure('visible','off'); end
    strain_figure(ff.epsyylon,res.XX,cfg.geom.N,cfg.geom.L,cfg.geom.thick,res.TT(end),f);
    if dosave
        fn = fullfile(opt.outdir, [cfg.name tag '_strain.png']);
        print(f,'-dpng','-r110',fn);
        fprintf('figure: %s\n', fn);
    end
    if ~show, close(f); else, drawnow; end
end
if ~isequal(opt.animate,false) && ~isempty(opt.animate)
    am = opt.animate;  if islogical(am) || isnumeric(am), am = ternary_(logical(am),'play','none'); end
    if ~any(strcmp(am,{'play','gif','both','none'})) && ~strcmp(am,'true')
        error('run_case: option ''animate'' must be false, true, ''gif'' or ''both''');
    end
    if strcmp(am,'true'), am = 'play'; end
    if ~strcmp(am,'none')
        ao = {'show', any(strcmp(am,{'play','both'}))};
        if any(strcmp(am,{'gif','both'}))
            if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
            ao = [ao {'gif', fullfile(opt.outdir,[cfg.name tag '_animation.gif'])}];
        end
        res.animation = animate_beam(res, cfg, ao{:}, opt.animate_opts{:});
        fprintf('animation done in %.1f s (%d frames)\n', res.animation.runtime_s, res.animation.nframes);
    end
end
if opt.save_mat
    if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
    TT = res.TT; XX = res.XX; midspan = res.midspan; midspan_node = res.midspan_node; x_midspan_node = res.x_midspan_node; peaks = res.peaks; %#ok<NASGU>
    save('-v7', fullfile(opt.outdir,[cfg.name tag '_result.mat']), 'TT','XX','midspan','midspan_node','x_midspan_node','peaks','cfg');
end
end

function v = ternary_(c,a,b)
if c, v = a; else, v = b; end
end

function v = getfield_or(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
