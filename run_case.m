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
%                   and the deformed beam + strain map (src/strain_figure.m, <case>_strain.png), and the martensite
%                   fractions versus time and position along the beam (src/martensite_xt.m, <case>_martensite_xt.png)
%     'dashboard'   true (default) | false: with 'figures' = 'show' / 'both' put the result plots in ONE window with tab
%                   buttons (src/dashboard.m: deflection, SMA vs linear, strain, martensite) instead of separate
%                   windows; false = the old separate windows. The PNG files are the same either way.
%     'animate'     false (default) | true = play the beam vibration in a window (src/animate_beam.m) |
%                   'gif' = write <outdir>/<case>_animation.gif (no window) | 'both' = window AND GIF.
%                   Slow motion, 25 fps, deflection exaggerated; works for every case.
%     'animate_compare' false (default) | true | 'gif' | 'both': the same for the SMA beam TOGETHER WITH the linear
%                   reference beam (src/animate_compare.m): beams overlaid + both mid-span histories with moving
%                   markers; GIF = <outdir>/<case>_animation_compare.gif. Runs the linear reference if needed.
%     'animate_fields' false (default) | true | 'gif' | 'both': animation of the martensite fractions xi_s+, xi_s-, xi_T over
%                   length x thickness while the beam rings down (src/animate_fields.m, ~100 snapshots recorded during
%                   the solve); GIF = <outdir>/<case>_fields_animation.gif. 'fields_opts': extra options (cell).
%     'animate_opts' cell of extra name/value options for animate_beam, e.g. {'speed',0.05,'strain',false}
%     'anim_speed'  playback speed of the animations, simulated seconds per real second (1 = real time, 0.1 = 10x slow
%                   motion; default 0.1): the start value of the speed control of the on-screen animations (adjustable
%                   between 0.02 and 1 while they play) and the speed of the GIF files. Default = unchanged GIF speed.
%     'anim_block'  true (default) | false. On-screen animations ('animate*' = true or 'both') are tabs of the results
%                   window (Beam animation, SMA vs linear animation, Martensite animation) and loop until that window
%                   is closed or Stop is pressed; run_case returns then (GNU Octave has no timer, so the playback loop
%                   runs inside run_case, as its last step, after all files are written). false = open the window with
%                   the animations paused and return at once; Play in the window starts the loop.
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
    'save_mat',true,'verbose',[],'plot_fields',false,'p',[],'animate',false,'animate_opts',{{}},'animate_compare',false,'animate_fields',false,'fields_opts',{{}},'dashboard',true,'anim_speed',[],'anim_block',true);
for k = 1:2:numel(varargin)
    if ~isfield(opt,varargin{k}), error('run_case: unknown option "%s"', varargin{k}); end
    opt.(varargin{k}) = varargin{k+1};
end
if ~isempty(opt.anim_speed) && ~(isnumeric(opt.anim_speed) && isscalar(opt.anim_speed) && opt.anim_speed > 0)
    error('run_case: option ''anim_speed'' must be a positive number (simulated seconds per real second, e.g. 0.1)');
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
want_fields_anim = ~isequal(opt.animate_fields,false) && ~isempty(opt.animate_fields);
if (opt.plot_fields || want_fields_anim) && ~strcmp(cfg.model,'linear') && ~(isfield(cfg.post,'n_snapshots') && cfg.post.n_snapshots > 0)
    cfg.post.n_snapshots = 100;                    % martensite fractions at ~100 times (src/sma_beam_solve.m): xt map, fields animation
end

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

ref = [];  refc = [];
if opt.compare_hle || ~isequal(opt.animate_compare,false)    % the comparison animation needs the linear reference too
    cr = case_hle_linear(); cr.time.Tf = cfg.time.Tf;
    refc = linear_beam_solve(cr);  refc.cfg = cr;
    res.reference = refc;
    if opt.compare_hle, ref = refc; end                       % ref = what the figures overlay
end

show = any(strcmp(opt.figures,{'show','both'}));
dosave = any(strcmp(opt.figures,{'save','both'}));
% ONE tabbed window (src/dashboard.m) instead of separate figure windows when figures are shown; the PNG files are
% written from hidden figures exactly as before. Falls back to separate windows if the toolkit has no uipanel.
use_dash = show && opt.dashboard && dashboard_supported();
vis_sep  = show && ~use_dash;                       % separate visible windows
if ~strcmp(opt.figures,'none') && (dosave || vis_sep)
    if vis_sep, f = figure; else, f = figure('visible','off'); end
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
        fprintf('figure: %s\n', fn);
    end
    if vis_sep, drawnow; else, close(f); end
end
if opt.plot_fields && isfield(res,'fields') && strcmp(cfg.model,'two_phase') && ~strcmp(opt.figures,'none') && (dosave || vis_sep)
    ff = res.fields;
    if ~exist(opt.outdir,'dir') && dosave, mkdir(opt.outdir); end
    % 1) martensite fractions over the beam at the last time step (thesis 'Simulazione_lunga')
    if vis_sep, f = figure; else, f = figure('visible','off'); end
    meshare(ff.epsyylon,ff.xysp,ff.xysm,ff.xyT,res.TT,cfg.geom.thick,cfg.geom.L,f,struct('ylabel',sma_section_label(cfg,'axis'),'suffix',sma_section_label(cfg,'suffix')));
    if dosave
        fn = fullfile(opt.outdir, [cfg.name tag '_fields.png']);
        print(f,'-dpng','-r110',fn);
        fprintf('figure: %s\n', fn);
    end
    if ~vis_sep, close(f); else, drawnow; end
    % 2) thesis-style deformed beam + strain map at the last time step
    if vis_sep, f = figure; else, f = figure('visible','off'); end
    strain_figure(ff.epsyylon,res.XX,cfg.geom.N,cfg.geom.L,cfg.geom.thick,res.TT(end),f,struct('ylabel',sma_section_label(cfg,'thickness'),'suffix',sma_section_label(cfg,'suffix')));
    if dosave
        fn = fullfile(opt.outdir, [cfg.name tag '_strain.png']);
        print(f,'-dpng','-r110',fn);
        fprintf('figure: %s\n', fn);
    end
    if ~vis_sep, close(f); else, drawnow; end
end
if opt.plot_fields && isfield(res,'snapshots') && ~strcmp(opt.figures,'none') && (dosave || vis_sep)
    % 3) martensite fractions versus time and position along the beam (thickness-averaged colour maps)
    if vis_sep, f = figure; else, f = figure('visible','off'); end
    set(f,'units','pixels','position',[100 80 700 760],'color','w','paperpositionmode','auto');
    martensite_xt(res, cfg, f);
    if dosave
        if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
        fn = fullfile(opt.outdir, [cfg.name tag '_martensite_xt.png']);
        print(f,'-dpng','-r110',fn);
        fprintf('figure: %s\n', fn);
    end
    if ~vis_sep, close(f); else, drawnow; end
end
[am_b,show_b,gif_b] = anim_mode_(opt.animate,'animate');
[am_c,show_c,gif_c] = anim_mode_(opt.animate_compare,'animate_compare');
[am_f,show_f,gif_f] = anim_mode_(opt.animate_fields,'animate_fields');
if show_f && ~isfield(res,'snapshots'), show_f = false; end          % linear beam: no martensite (message below)
sp_screen = anim_start_speed_(opt);                                    % start speed of the on-screen animations
% on-screen animations are tabs of the dashboard (one window); without uipanel support they open their own windows
dash_anim = opt.dashboard && dashboard_supported() && (show_b || show_c || show_f);
res.dashboard = [];
if use_dash || dash_anim
    try
        A = struct('beam',show_b,'compare',show_c,'fields',show_f,'refc',refc,'speed',sp_screen, ...
            'beam_opts',{anim_opts_(opt.animate_opts)},'fields_opts',{anim_opts_(opt.fields_opts)},'autoplay',opt.anim_block);
        res.dashboard = dashboard(res, cfg, ref, 'title', [cfg.name tag], 'plots', use_dash, 'anim', A);
        drawnow;
    catch err
        warning('run_case: the dashboard window could not be created (%s); use ''dashboard'',false for separate windows', err.message);
        dash_anim = false;  res.dashboard = [];
    end
end
if isempty(res.dashboard), res = rmfield(res,'dashboard'); end
if isfield(res,'dashboard') && dash_anim
    show_b = false;  show_c = false;  show_f = false;                   % already in the dashboard: nothing else is opened
end
gs = {};  if ~isempty(opt.anim_speed), gs = {'speed',opt.anim_speed}; end    % GIF speed
if ~strcmp(am_b,'none')
    ao = [{'show', show_b}, gs];
    if gif_b
        if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
        ao = [ao {'gif', fullfile(opt.outdir,[cfg.name tag '_animation.gif'])}];
    end
    if gif_b || show_b
        res.animation = animate_beam(res, cfg, ao{:}, opt.animate_opts{:});
        fprintf('animation done in %.1f s (%d frames)\n', res.animation.runtime_s, res.animation.nframes);
    end
end
if ~strcmp(am_c,'none')
    ao = [{'show', show_c}, gs];
    if gif_c
        if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
        ao = [ao {'gif', fullfile(opt.outdir,[cfg.name tag '_animation_compare.gif'])}];
    end
    if gif_c || show_c
        res.animation_compare = animate_compare(res, refc, cfg, refc.cfg, ao{:}, opt.animate_opts{:});
        fprintf('comparison animation done in %.1f s (%d frames)\n', res.animation_compare.runtime_s, res.animation_compare.nframes);
    end
end
if ~strcmp(am_f,'none')
    if ~isfield(res,'snapshots')
        fprintf('animate_fields: no martensite fractions in this model (linear beam) - skipped\n');
    elseif gif_f || show_f
        ao = [{'show', show_f}, gs];
        if gif_f
            if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
            ao = [ao {'gif', fullfile(opt.outdir,[cfg.name tag '_fields_animation.gif'])}];
        end
        res.animation_fields = animate_fields(res, cfg, ao{:}, opt.fields_opts{:});
        fprintf('fields animation done in %.1f s (%d frames)\n', res.animation_fields.runtime_s, res.animation_fields.nframes);
    end
end
if opt.save_mat
    if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
    TT = res.TT; XX = res.XX; midspan = res.midspan; midspan_node = res.midspan_node; x_midspan_node = res.x_midspan_node; peaks = res.peaks; %#ok<NASGU>
    save('-v7', fullfile(opt.outdir,[cfg.name tag '_result.mat']), 'TT','XX','midspan','midspan_node','x_midspan_node','peaks','cfg');
end
% on-screen animations: loop until the window is closed or Stop is pressed (everything is saved at this point)
if dash_anim && isfield(res,'dashboard') && ishandle(res.dashboard) && opt.anim_block
    fprintf('animations are playing in the results window - close it or press Stop to continue\n');
    dashboard_play(res.dashboard);
end
end

function sp = anim_start_speed_(opt)
% start speed of the on-screen animations: 'anim_speed', else a numeric 'speed' in animate_opts, else 0.1
sp = 0.1;
if ~isempty(opt.anim_speed), sp = opt.anim_speed; return; end
c = opt.animate_opts;
for k = 1:2:numel(c)-1
    if ischar(c{k}) && strcmp(c{k},'speed') && isnumeric(c{k+1}) && isscalar(c{k+1}), sp = c{k+1}; end
end
end

function c = anim_opts_(c)
% extra options handed to the dashboard animation tabs: drop what only applies to files / own windows
keep = true(size(c));
for k = 1:2:numel(c)-1
    if ischar(c{k}) && any(strcmp(c{k},{'gif','show','parent','still','speed'})), keep(k:k+1) = false; end
end
c = c(keep);
end

function [mode,doshow,dogif] = anim_mode_(v,name)
% value of an 'animate*' option -> 'none' | 'play' | 'gif' | 'both'
if isempty(v) || isequal(v,false) || isequal(v,0), mode = 'none';
elseif islogical(v) || isnumeric(v), mode = 'play';
elseif strcmp(v,'true'), mode = 'play';
elseif any(strcmp(v,{'play','gif','both','none'})), mode = v;
else, error('run_case: option ''%s'' must be false, true, ''gif'' or ''both''', name);
end
doshow = any(strcmp(mode,{'play','both'}));  dogif = any(strcmp(mode,{'gif','both'}));
end

function v = getfield_or(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
