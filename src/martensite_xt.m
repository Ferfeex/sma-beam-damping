function martensite_xt(res, cfg, parent, varargin)
% MARTENSITE_XT  evolution of the martensite fractions in time and along the beam, as colour maps.
%
%   martensite_xt(res, cfg, parent, 'name', value, ...)
%     res     result of run_case with snapshots (cfg.post.n_snapshots > 0; run_case('...','plot_fields',true))
%     parent  figure (or uitab) to draw in
%   One map per fraction, horizontal axis = position x along the beam, vertical axis = time t (running upwards):
%     xi_s+  martensite induced by tension,   xi_s-  induced by compression,   xi_T  twinned martensite
%   (two_phase models; the super-elastic pre-strain model has only xi_s+). Every map has its own colour scale.
%   Options:
%     'fibre'  'mean' (default; average over the section thickness) | 'top' | 'bottom' | 'max'
%   For the sandwich beam the maps refer to the SMA core only (the GFRP skins are linear and have no martensite).
%   The beam is released from the first mid-span maximum, so t = 0 is the first maximum of the mid-span deflection.
o = struct('fibre','mean');
for k = 1:2:numel(varargin)
    if ~isfield(o,varargin{k}), error('martensite_xt: unknown option "%s"', varargin{k}); end
    o.(varargin{k}) = varargin{k+1};
end
if ~isfield(res,'snapshots'), error('martensite_xt: the result has no snapshots (set cfg.post.n_snapshots)'); end
names = {'xisp'};  labs = {'\xi_s^+  (tension)'};
if ~isempty(res.snapshots.xism)
    names = {'xisp','xism','xiT'};  labs = {'\xi_s^+  (tension)','\xi_s^-  (compression)','\xi_T  (twinned)'};
end
fibres = struct('mean','section average','top','top fibre (y = +h/2)','bottom','bottom fibre (y = -h/2)','max','section maximum');
n = numel(names);
top = 0.90;  bot = 0.10;  gap = 0.06;  h = (top-bot-gap*(n-1))/n;
for i = 1:n
    [M,x,t] = snapshot_fields(res, cfg, names{i}, o.fibre);
    yb = top - i*h - (i-1)*gap;
    ax = axes('parent',parent,'position',[0.11 yb 0.72 h]);
    if numel(t) > 1                              % snapshot times are only nearly even: resample on an even grid
        tt = linspace(t(1),t(end),numel(t)).';  M = interp1(t,M,tt,'linear');
    else, tt = [t t+eps]; M = [M; M]; end
    imagesc(ax,x, tt, M); set(ax,'ydir','normal');
    lo = min(M(:));  hi = max(M(:));  if hi-lo < 1e-12, hi = lo + 1e-12; end
    colormap(ax,jet(128)); set(ax,'clim',[lo hi]);
    ylabel(ax,'time t [s]');
    if i < n, set(ax,'xticklabel',[]); else, xlabel(ax,'position along the beam x [mm]'); end
    title(ax, sprintf('%s - %s%s', labs{i}, fibres.(o.fibre), sma_section_label(cfg,'suffix')), 'fontsize',10);
    hand_colorbar(parent,[0.86 yb 0.025 h],[lo hi],jet(128),'',  '%.4f');
end
end
