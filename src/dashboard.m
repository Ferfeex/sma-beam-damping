function fig = dashboard(res, cfg, ref, varargin)
% DASHBOARD  ONE window with a row of tab buttons that gathers the result plots of a run
%   (instead of a pile of overlapping figure windows):
%     Deflection         mid-span displacement versus time (with the linear reference if it was run)
%     SMA vs linear      (only with the linear reference) both mid-span curves and the decay of their peaks
%     Strain             deformed beam and strain map at the last time step         (plot_fields, two-phase models)
%     Martensite (final) xi_s+, xi_s-, xi_T over length x thickness at the last step (plot_fields, two-phase models)
%     Martensite (t, x)  xi versus time and position along the beam                 (plot_fields)
%     Beam animation, SMA vs linear animation, Martensite animation
%                        looping on-screen animations with Play/Pause, speed slider, -/+ and Stop controls
%                        (only the visible tab is animated; played by dashboard_play.m)
%
%   fig = dashboard(res, cfg, ref, 'name', value, ...)
%     res, cfg   result and configuration of run_case (res.fields / res.snapshots exist after 'plot_fields',true)
%     ref        result of the linear reference run ([] = none)
%     'title'    window title suffix (default: case name)
%     'visible'  true (default) | false (build it without showing it, used by the tests)
%     'plots'    true (default) | false: include the result-plot tabs (false = animation tabs only)
%     'anim'     struct with the animation tabs to build: .beam, .compare, .fields (logical), .refc (linear reference
%                result for the comparison tab), .speed (start speed, simulated s per real s, 0.02..1),
%                .beam_opts / .fields_opts (cells of extra animate_beam / animate_fields options), .autoplay
%                (true: the animations start playing; false: they wait for Play)
%
%   The tabs are ordinary uicontrol buttons + uipanels, because GNU Octave (10.3) has no uitabgroup; this works
%   unchanged in the Octave Qt and FLTK toolkits and in MATLAB. Needs a graphics toolkit with uipanel support
%   (dashboard_supported). The window can be resized; the animations stay in their own windows.
o = struct('title','','visible',true,'plots',true,'anim',struct());
for k = 1:2:numel(varargin)
    if ~isfield(o,varargin{k}), error('dashboard: unknown option "%s"', varargin{k}); end
    o.(varargin{k}) = varargin{k+1};
end
if isempty(o.title), o.title = strrep(cfg.name,'_','\_'); end
if nargin < 3, ref = []; end
A = struct('beam',false,'compare',false,'fields',false,'refc',[],'speed',0.1,'beam_opts',{{}},'fields_opts',{{}},'autoplay',true);
fa = fieldnames(o.anim);
for k = 1:numel(fa), A.(fa{k}) = o.anim.(fa{k}); end

names = {};  draw = {};  isanim = [];
if o.plots
    names{end+1} = 'Deflection';          draw{end+1} = @(p) tab_deflection(p, res, cfg, ref);  isanim(end+1) = false;
    if ~isempty(ref)
        names{end+1} = 'SMA vs linear';   draw{end+1} = @(p) tab_compare(p, res, ref, cfg);     isanim(end+1) = false;
    end
    if isfield(res,'fields') && strcmp(cfg.model,'two_phase')
        names{end+1} = 'Strain';              draw{end+1} = @(p) strain_figure(res.fields.epsyylon,res.XX,cfg.geom.N,cfg.geom.L,cfg.geom.thick,res.TT(end),p,struct('ylabel',sma_section_label(cfg,'thickness'),'suffix',sma_section_label(cfg,'suffix')));  isanim(end+1) = false;
        names{end+1} = 'Martensite (final)';  draw{end+1} = @(p) tab_fields_final(p, res, cfg);  isanim(end+1) = false;
    end
    if isfield(res,'snapshots')
        names{end+1} = 'Martensite (t, x)';   draw{end+1} = @(p) martensite_xt(res, cfg, p);     isanim(end+1) = false;
    end
end
sp0 = min(1, max(0.02, A.speed));
if A.beam
    names{end+1} = 'Beam animation';
    draw{end+1} = @(p) animate_beam(res, cfg, 'parent',p, 'verbose',false, 'speed',sp0, A.beam_opts{:});  isanim(end+1) = true;
end
if A.compare && ~isempty(A.refc)
    names{end+1} = 'SMA vs linear animation';
    draw{end+1} = @(p) animate_compare(res, A.refc, cfg, A.refc.cfg, 'parent',p, 'verbose',false, 'speed',sp0);  isanim(end+1) = true;
end
if A.fields && isfield(res,'snapshots')
    names{end+1} = 'Martensite animation';
    draw{end+1} = @(p) animate_fields(res, cfg, 'parent',p, 'verbose',false, A.fields_opts{:});  isanim(end+1) = true;
end
if isempty(names), error('dashboard: nothing to show'); end

vis = 'off'; if o.visible, vis = 'on'; end
fig = figure('visible',vis,'name',['SMA beam - ' strrep(o.title,'\_','_') ' - results'],'numbertitle','off', ...
    'units','pixels','position',[80 60 940 760],'color',[1 1 1],'tag','sma_dashboard');
n = numel(names);
% tab buttons: one row; plots and animations in two rows when there are more than 5 tabs
if n > 5 && any(isanim) && any(~isanim), rowof = 1 + double(isanim); else, rowof = ones(1,n); end
nrow = max(rowof);  rh = 0.052;  ptop = 1 - nrow*(rh+0.006) - 0.004;
panels = zeros(1,n);  btns = zeros(1,n);
for r = 1:nrow
    idx = find(rowof == r);  bw = min(0.22, 0.96/numel(idx));
    for q = 1:numel(idx)
        k = idx(q);
        btns(k) = uicontrol(fig,'style','togglebutton','units','normalized','position',[0.005+(q-1)*(bw+0.005) 1-r*(rh+0.006) bw rh], ...
            'string',names{k},'fontweight','bold','fontsize',9,'value',0,'callback',@(h,e) dashboard_select(fig,k));
    end
end
for k = 1:n
    panels(k) = uipanel('parent',fig,'units','normalized','position',[0 0 1 ptop],'bordertype','none', ...
        'backgroundcolor',[1 1 1],'visible','off');
end
players = cell(1,n);  ctrl = zeros(n,4);
for k = 1:n
    if isanim(k)
        [canvas, ctrl(k,:)] = anim_controls(fig, panels(k), k, sp0, A.autoplay);
        out = draw{k}(canvas);
        players{k} = out.player;
    else
        draw{k}(panels(k));
    end
end
setappdata(fig,'dash',struct('panels',panels,'buttons',btns,'names',{names},'isanim',isanim));
setappdata(fig,'play_players',players);
setappdata(fig,'play_ctrl',ctrl);                       % [play button, slider, speed text, stop button] per tab
setappdata(fig,'play_speed',sp0*ones(1,n));
setappdata(fig,'play_on',logical(A.autoplay)*ones(1,n));
setappdata(fig,'play_run',false);
setappdata(fig,'play_phase',zeros(1,n));
first = find(isanim,1);  if isempty(first), first = 1; end   % the first animation tab is shown first
dashboard_select(fig,first);
if o.visible, drawnow; end
end

% ======================================================================================
function [canvas, h] = anim_controls(fig, p, k, sp0, autoplay)
% canvas panel for the animation plus a control strip: Play/Pause, -, slider, +, speed text, Stop
canvas = uipanel('parent',p,'units','normalized','position',[0 0.075 1 0.925],'bordertype','none','backgroundcolor',[1 1 1]);
bar = uipanel('parent',p,'units','normalized','position',[0 0 1 0.075],'bordertype','etchedin','backgroundcolor',[0.94 0.94 0.94]);
lbl = 'Pause'; if ~autoplay, lbl = 'Play'; end
h = zeros(1,4);
h(1) = uicontrol(bar,'style','pushbutton','units','normalized','position',[0.01 0.12 0.10 0.76],'string',lbl,'fontweight','bold', ...
    'tooltipstring','play / pause this animation','callback',@(s,e) dashboard_ctrl(fig,k,'toggle'));
uicontrol(bar,'style','text','units','normalized','position',[0.13 0.15 0.08 0.6],'string','speed','backgroundcolor',[0.94 0.94 0.94],'horizontalalignment','right');
uicontrol(bar,'style','pushbutton','units','normalized','position',[0.22 0.12 0.04 0.76],'string','-','fontweight','bold', ...
    'tooltipstring','slower','callback',@(s,e) dashboard_ctrl(fig,k,'slower'));
h(2) = uicontrol(bar,'style','slider','units','normalized','position',[0.265 0.2 0.30 0.6],'min',log10(0.02),'max',0,'value',log10(sp0), ...
    'tooltipstring','playback speed, 0.02x - 1x real time','callback',@(s,e) dashboard_ctrl(fig,k,'slider'));
uicontrol(bar,'style','pushbutton','units','normalized','position',[0.57 0.12 0.04 0.76],'string','+','fontweight','bold', ...
    'tooltipstring','faster','callback',@(s,e) dashboard_ctrl(fig,k,'faster'));
h(3) = uicontrol(bar,'style','text','units','normalized','fontsize',9,'position',[0.615 0.15 0.25 0.6],'string',dashboard_ctrl(fig,k,'label',sp0), ...
    'backgroundcolor',[0.94 0.94 0.94],'horizontalalignment','left');
h(4) = uicontrol(bar,'style','pushbutton','units','normalized','position',[0.87 0.12 0.12 0.76],'string','Stop','fontweight','bold', ...
    'tooltipstring','stop the playback loop (the window stays open; Play restarts it)','callback',@(s,e) dashboard_ctrl(fig,k,'stop'));
end

% ======================================================================================
function tab_deflection(p, res, cfg, ref)
ax = axes('parent',p,'position',[0.08 0.09 0.89 0.83]);
plot(ax, res.TT, res.midspan, 'LineWidth', 1.5); hold(ax,'on')
if ~isempty(ref), plot(ax, ref.TT, ref.midspan, 'LineWidth', 1.5); end
grid(ax,'on'); xlabel(ax,'Time t [s]'); ylabel(ax,'Displacement \delta [mm]');
ttl = sprintf('Mid-span vertical displacement (x = %g mm) - %s', cfg.geom.L/2, strrep(cfg.name,'_','\_'));
if isfield(cfg,'p_modified') && cfg.p_modified, ttl = [ttl ' (modified p)']; end
title(ax,ttl);
if ~isempty(ref)
    legend(ax,sprintf('SMA (loss-factor damping), x = %g mm',cfg.geom.L/2),sprintf('linear reference, E = 26 GPa, Rayleigh 1.3 %%, x = %g mm',cfg.geom.L/2));
end
end

function tab_compare(p, res, ref, cfg)
c1 = [0 0.2 0.8];  c2 = [0.85 0.1 0.1];
ax1 = axes('parent',p,'position',[0.08 0.56 0.89 0.37]);
plot(ax1, ref.TT, ref.midspan, '-', 'color',c2, 'linewidth',1.2); hold(ax1,'on')
plot(ax1, res.TT, res.midspan, '-', 'color',c1, 'linewidth',1.2);
grid(ax1,'on'); ylabel(ax1,sprintf('w(x=%g mm) [mm]',cfg.geom.L/2)); xlabel(ax1,'time t [s]');
title(ax1,'SMA beam against the linear reference beam');
legend(ax1,'linear reference (Rayleigh)','SMA (loss-factor damping)','location','northeast');
[t1,x1] = sma_peaks(res.TT,res.midspan);  [t2,x2] = sma_peaks(ref.TT,ref.midspan);
ax2 = axes('parent',p,'position',[0.08 0.08 0.89 0.37]);
semilogy(ax2, t2, x2, 's-', 'color',c2, 'linewidth',1.2, 'markersize',5, 'markerfacecolor',c2); hold(ax2,'on')
semilogy(ax2, t1, x1, 'o-', 'color',c1, 'linewidth',1.2, 'markersize',5, 'markerfacecolor',c1);
grid(ax2,'on'); xlabel(ax2,'time t [s]'); ylabel(ax2,'height of the mid-span peaks [mm]');
title(ax2,'Decay of the oscillation (peaks, logarithmic scale)');
end

function tab_fields_final(p, res, cfg)
% xi_s+, xi_s-, xi_T over length x thickness at the last time step (same data and smoothing as src/meshare.m)
ff = res.fields;  thick = cfg.geom.thick;  L = cfg.geom.L;
y = linspace(thick/2,-thick/2,10)';
vals = {ff.xysp, ff.xysm, ff.xyT};  labs = {'\xi_s^+ (tension)','\xi_s^- (compression)','\xi_T (twinned)'};
x = linspace(0,L,size(reshape(vals{1}(:,end),numel(y),[]),2))';
[X,Y] = meshgrid(x,y);  [xx,yy] = meshgrid(linspace(0,L,400), linspace(thick/2,-thick/2,200));
top = 0.91;  bot = 0.09;  gap = 0.07;  h = (top-bot-2*gap)/3;
for i = 1:3
    F = interp2(X,Y,smooth_mean2d(reshape(vals{i}(:,end),numel(y),[]),1,8),xx,yy);
    yb = top - i*h - (i-1)*gap;
    ax = axes('parent',p,'position',[0.10 yb 0.72 h]);
    imagesc(ax, xx(1,:), yy(:,1), F); set(ax,'ydir','normal'); colormap(ax,jet(128));
    lo = min(F(:)); hi = max(F(:)); if hi-lo < 1e-12, hi = lo+1e-12; end
    set(ax,'clim',[lo hi]);
    ylabel(ax,sma_section_label(cfg,'axis')); title(ax,labs{i},'fontsize',10);
    if i < 3, set(ax,'xticklabel',[]); else, xlabel(ax,'length [mm]'); end
    if i == 1, title(ax,sprintf('t = %.3g s    %s%s', res.TT(end), labs{i}, sma_section_label(cfg,'suffix')),'fontsize',10); end
    hand_colorbar(p,[0.85 yb 0.025 h],[lo hi],jet(128),'','%.4f');
end
end
