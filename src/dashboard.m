function fig = dashboard(res, cfg, ref, varargin)
% DASHBOARD  ONE window with a row of tab buttons that gathers the result plots of a run
%   (instead of a pile of overlapping figure windows):
%     Deflection         mid-span displacement versus time (with the linear reference if it was run)
%     SMA vs linear      (only with the linear reference) both mid-span curves and the decay of their peaks
%     Strain             deformed beam and strain map at the last time step         (plot_fields, two-phase models)
%     Martensite (final) xi_s+, xi_s-, xi_T over length x thickness at the last step (plot_fields, two-phase models)
%     Martensite (t, x)  xi versus time and position along the beam                 (plot_fields)
%
%   fig = dashboard(res, cfg, ref, 'name', value, ...)
%     res, cfg   result and configuration of run_case (res.fields / res.snapshots exist after 'plot_fields',true)
%     ref        result of the linear reference run ([] = none)
%     'title'    window title suffix (default: case name)
%     'visible'  true (default) | false (build it without showing it, used by the tests)
%
%   The tabs are ordinary uicontrol buttons + uipanels, because GNU Octave (10.3) has no uitabgroup; this works
%   unchanged in the Octave Qt and FLTK toolkits and in MATLAB. Needs a graphics toolkit with uipanel support
%   (dashboard_supported). The window can be resized; the animations stay in their own windows.
o = struct('title','','visible',true);
for k = 1:2:numel(varargin)
    if ~isfield(o,varargin{k}), error('dashboard: unknown option "%s"', varargin{k}); end
    o.(varargin{k}) = varargin{k+1};
end
if isempty(o.title), o.title = strrep(cfg.name,'_','\_'); end
if nargin < 3, ref = []; end

names = {};  draw = {};
names{end+1} = 'Deflection';          draw{end+1} = @(p) tab_deflection(p, res, cfg, ref);
if ~isempty(ref)
    names{end+1} = 'SMA vs linear';   draw{end+1} = @(p) tab_compare(p, res, ref, cfg);
end
if isfield(res,'fields') && strcmp(cfg.model,'two_phase')
    names{end+1} = 'Strain';              draw{end+1} = @(p) strain_figure(res.fields.epsyylon,res.XX,cfg.geom.N,cfg.geom.L,cfg.geom.thick,res.TT(end),p);
    names{end+1} = 'Martensite (final)';  draw{end+1} = @(p) tab_fields_final(p, res, cfg);
end
if isfield(res,'snapshots')
    names{end+1} = 'Martensite (t, x)';   draw{end+1} = @(p) martensite_xt(res, cfg, p);
end

vis = 'off'; if o.visible, vis = 'on'; end
fig = figure('visible',vis,'name',['SMA beam - ' strrep(o.title,'\_','_') ' - results'],'numbertitle','off', ...
    'units','pixels','position',[80 60 940 760],'color',[1 1 1],'tag','sma_dashboard');
n = numel(names);  bw = min(0.22, 0.96/n);
panels = zeros(1,n);  btns = zeros(1,n);
for k = 1:n
    panels(k) = uipanel('parent',fig,'units','normalized','position',[0 0 1 0.93],'bordertype','none', ...
        'backgroundcolor',[1 1 1],'visible','off');
    btns(k) = uicontrol(fig,'style','togglebutton','units','normalized','position',[0.005+(k-1)*(bw+0.005) 0.94 bw 0.052], ...
        'string',names{k},'fontweight','bold','value',0,'callback',@(h,e) dashboard_select(fig,k));
end
for k = 1:n
    draw{k}(panels(k));
end
setappdata(fig,'dash',struct('panels',panels,'buttons',btns,'names',{names}));
dashboard_select(fig,1);
if o.visible, drawnow; end
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
    ylabel(ax,'h [mm]'); title(ax,labs{i},'fontsize',10);
    if i < 3, set(ax,'xticklabel',[]); else, xlabel(ax,'length [mm]'); end
    if i == 1, title(ax,sprintf('t = %.3g s    %s', res.TT(end), labs{i}),'fontsize',10); end
    hand_colorbar(p,[0.85 yb 0.025 h],[lo hi],jet(128),'','%.4f');
end
end
