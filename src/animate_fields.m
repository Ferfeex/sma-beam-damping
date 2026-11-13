function out = animate_fields(res, cfg, varargin)
% ANIMATE_FIELDS  animation of the martensite fractions xi_s+, xi_s-, xi_T over the beam (length x thickness)
%   as they build up during the vibration (and GIF export).
%
%   out = animate_fields(res, cfg, 'name', value, ...)
%     res must hold snapshots (cfg.post.n_snapshots > 0; run_case does it for 'animate_fields' / 'plot_fields').
%   The maps are the ones of the final-time figure (<case>_fields.png: same smoothing, src/smooth_mean2d.m),
%   with FIXED colour limits over the whole run (one per fraction), a time label and, at the bottom, the
%   mid-span deflection with a moving marker. The super-elastic pre-strain model shows xi_s+ only.
%   Options (defaults)
%     'gif'     ''     GIF file to write ('' = none); it loops forever
%     'show'    false  play in a visible window
%     'fps'     10     frames per second (one frame per snapshot, ~100 frames)
%     'size'    [560 150+150*maps] pixels
%     'step'    1      use every step-th snapshot only (smaller GIF)
%     'still'   ''     PNG with the last frame
%     'verbose' true
%   Returns out: .gif, .nframes, .t, .fps, .frame_size, .runtime_s, .clim
t_start = tic;
o = struct('gif','','show',false,'fps',10,'size',[],'step',1,'still','','verbose',true);
for k = 1:2:numel(varargin)
    if ~isfield(o,varargin{k}), error('animate_fields: unknown option "%s"', varargin{k}); end
    o.(varargin{k}) = varargin{k+1};
end
if ~isfield(res,'snapshots'), error('animate_fields: the result has no snapshots (set cfg.post.n_snapshots)'); end
S = res.snapshots;  L = cfg.geom.L;  thick = cfg.geom.thick;
names = {'xisp'};  labs = {'\xi_s^+ (tension)'};
if ~isempty(S.xism), names = {'xisp','xism','xiT'};  labs = {'\xi_s^+ (tension)','\xi_s^- (compression)','\xi_T (twinned)'}; end
nm = numel(names);  sel = 1:o.step:numel(S.t);  nfr = numel(sel);  t = S.t(sel).';  t = t(:);
if isempty(o.size), o.size = [560 150+150*nm]; end
ny = S.n_y;  nc = size(S.xisp,1)/ny;
ys_ = linspace(thick/2,-thick/2,ny).';  xs_ = linspace(0,L,nc);
[Xs,Ys] = meshgrid(xs_,ys_);
[xi,yi] = meshgrid(linspace(0,L,240), linspace(thick/2,-thick/2,60));
G = cell(nm,1);  clim = zeros(nm,2);
for i = 1:nm                                          % smoothed, interpolated maps of every frame (60 x 240 x nfr)
    G{i} = zeros(size(xi,1), size(xi,2), nfr);
    for j = 1:nfr
        G{i}(:,:,j) = interp2(Xs,Ys,smooth_mean2d(reshape(S.(names{i})(:,sel(j)),ny,nc),1,8),xi,yi);
    end
    clim(i,:) = [min(G{i}(:)) max(G{i}(:))];
    if diff(clim(i,:)) < 1e-12, clim(i,2) = clim(i,1) + 1e-12; end
end
[Tu,iu] = unique(res.TT(:));  mid = interp1(Tu, res.midspan(iu), t, 'linear', 'extrap');

if o.show, fig = figure('name','martensite evolution','numbertitle','off'); else, fig = figure('visible','off'); end
set(fig,'units','pixels','position',[100 100 o.size],'color','w','paperpositionmode','auto');
np = nm + 1;  gap = 0.045;  top = 0.94;  bot = 0.08;  hp = (top-bot-gap*(np-1))/np;
him = zeros(nm,1);
for i = 1:nm
    yb = top - i*hp - (i-1)*gap;
    ax = axes('position',[0.11 yb 0.70 hp]);
    him(i) = imagesc(xi(1,:), yi(:,1), G{i}(:,:,1)); set(ax,'ydir','normal'); colormap(ax,jet(128)); set(ax,'clim',clim(i,:));
    ylabel('h [mm]');  title(labs{i},'fontsize',9);
    set(ax,'xticklabel',[]);
    hand_colorbar(fig,[0.84 yb 0.025 hp],clim(i,:),jet(128),'','%.4f');
    if i == 1, htt = get(ax,'title'); end
end
yb = bot;
ax3 = axes('position',[0.11 yb 0.70 hp]); hold on
plot(Tu, res.midspan(iu), 'color',[0.6 0.6 0.6]);
hm = plot(t(1),mid(1),'ro','markerfacecolor','r','markersize',6);
axis([t(1) t(end) 1.1*min(res.midspan) 1.1*max(res.midspan)]); grid on; box on
xlabel('time t [s]'); ylabel(sprintf('w(x=%g mm) [mm]',L/2));
upd = @(j) local_update(j, him, G, hm, t, mid, htt, labs{1});

dogif = ~isempty(o.gif);  frames = cell(nfr,1);
if dogif
    sidx = unique(round(linspace(1,nfr,min(nfr,7))));  Sx = [];
    for j = sidx, upd(j); Sx = [Sx; reshape(gif_grab_frame(fig),[],3)]; end %#ok<AGROW>
    [map,lut] = gif_palette(Sx);
end
for j = 1:nfr
    upd(j);
    if dogif, frames{j} = gif_quantize(gif_grab_frame(fig), lut); elseif o.show, pause(1/o.fps); end
    if o.verbose && mod(j,50)==0, fprintf('  fields frame %d / %d\n', j, nfr); end
end
if ~isempty(o.still), upd(nfr); imwrite(gif_grab_frame(fig), o.still); end
if o.show, drawnow; else, close(fig); end
out = struct('gif','','nframes',nfr,'t',t,'fps',o.fps,'frame_size',[],'runtime_s',0,'clim',clim);
if dogif
    [h,w] = size(frames{1});
    gif_write(o.gif, frames, map, o.fps);
    out.gif = o.gif;  out.frame_size = [w h];
    if o.verbose
        fi = dir(o.gif);
        fprintf('fields animation: %d frames, %.1f MB -> %s\n', nfr, fi.bytes/1e6, o.gif);
    end
end
out.runtime_s = toc(t_start);
end

function local_update(j, him, G, hm, t, mid, htt, lab1)
for i = 1:numel(him), set(him(i),'cdata',G{i}(:,:,j)); end
set(hm,'xdata',t(j),'ydata',mid(j));
set(htt,'string',sprintf('t = %.3f s    %s',t(j),lab1),'fontsize',9);
drawnow;
end
