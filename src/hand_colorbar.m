function ax = hand_colorbar(parent, pos, clim, cmap, label, fmt, scale)
% HAND_COLORBAR  small colour bar drawn as an ordinary axes at the normalized position pos = [x y w h] of parent
%   (figure or uitab). Same look in Octave (any toolkit) and MATLAB, and it works inside tabs, where 'colorbar'
%   is not reliable. Tick labels are v*scale printed with fmt (default '%.3g', scale 1).
if nargin < 6 || isempty(fmt), fmt = '%.3g'; end
if nargin < 7, scale = 1; end
ax = axes('parent',parent,'position',pos);
v = linspace(clim(1),clim(2),128);
imagesc(ax,[0 1], v, v(:)); set(ax,'ydir','normal');
colormap(ax,cmap); set(ax,'clim',clim);
tk = linspace(clim(1),clim(2),5);
set(ax,'xtick',[],'yaxislocation','right','fontsize',9,'ytick',tk, ...
    'yticklabel',arrayfun(@(q) sprintf(fmt,q*scale),tk,'uniformoutput',false));
if nargin >= 5 && ~isempty(label), title(ax,label,'fontsize',9); end
end
