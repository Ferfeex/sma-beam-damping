function out = animate_compare(res, ref, cfg, cfg_ref, varargin)
% ANIMATE_COMPARE  animation of the SMA beam and the linear reference beam vibrating TOGETHER (and GIF export).
%
%   out = animate_compare(res, ref, cfg, cfg_ref, 'name', value, ...)
%       res, cfg      result / configuration of the SMA run (run_case)
%       ref, cfg_ref  result / configuration of the linear reference (case_hle_linear, run by run_case when
%                     'compare_hle' or 'animate_compare' is used)
%
%   Two panels with the SAME time axis and the same deflection scale:
%     top     both beams overlaid: SMA (blue), linear reference (red), deflection exaggerated (factor in the label)
%     bottom  mid-span deflection (x = L/2) of both versus time, with a moving marker on each curve
%   The SMA beam dies out faster because every stress reversal dissipates energy; the linear beam only has
%   the small Rayleigh damping.
%
%   Options (defaults)
%     'gif'          ''     GIF file to write ('' = none); it loops forever (src/gif_write.m)
%     'show'         false  play in a visible window
%     'fps'          25
%     'speed'        'auto' simulated seconds per real second (auto = 0.1: 10x slow motion)
%     'exaggeration' 'auto'
%     'tspan'        []     [t0 t1] (default: the common time span of both runs)
%     'maxframes'    600
%     'size'         [560 480]
%     'still'        ''     PNG with one still frame (at the middle of the first quarter of the run)
%     'verbose'      true
%     'parent'       []     figure or uipanel to build the animation in (dashboard tab); nothing is played or written, out.player
%                           holds the frame update function (see dashboard_play.m)
%     'playframes'   300    number of uniform frame times prepared for the on-screen player
%   Returns out: .gif, .nframes, .t, .exaggeration, .speed, .fps, .frame_size, .runtime_s
t_start = tic;
o = struct('gif','','show',false,'fps',25,'speed','auto','exaggeration','auto','tspan',[],'maxframes',600, ...
           'size',[560 480],'still','','verbose',true,'strain',false,'history',true,'parent',[],'playframes',300);   % strain/history: accepted and ignored, so that one 'animate_opts' list can serve animate_beam too
for k = 1:2:numel(varargin)
    if ~isfield(o,varargin{k}), error('animate_compare: unknown option "%s"', varargin{k}); end
    o.(varargin{k}) = varargin{k+1};
end
% ---- common frame times -------------------------------------------------------------------------
[T1,i1] = unique(res.TT(:));  [T2,i2] = unique(ref.TT(:));
t0 = max([T1(1) T2(1)]);  t1 = min([T1(end) T2(end)]);
if ~isempty(o.tspan), t0 = max(t0,o.tspan(1)); t1 = min(t1,o.tspan(2)); end
if t1 <= t0, error('animate_compare: empty time span'); end
if ischar(o.speed), speed = 0.1; else, speed = o.speed; end
inparent = ~isempty(o.parent);
nfr = max(2, floor((t1-t0)*o.fps/speed) + 1);
if inparent, nfr = o.playframes; end
if ~inparent && nfr > o.maxframes
    if o.verbose, fprintf('animate_compare: %d frames needed, capped at %d\n', nfr, o.maxframes); end
    nfr = o.maxframes;  speed = (t1-t0)*o.fps/(nfr-1);
end
tf = linspace(t0,t1,nfr).';
N1 = cfg.geom.N;  N2 = cfg_ref.geom.N;  L = cfg.geom.L;
U1 = interp1(T1, res.XX(i1,2*N1-1:end), tf, 'linear');
U2 = interp1(T2, ref.XX(i2,2*N2-1:end), tf, 'linear');
m1 = interp1(T1, res.midspan(i1), tf, 'linear');
m2 = interp1(T2, ref.midspan(i2), tf, 'linear');
[xg,H1] = beam_shape_matrix(N1, L, 201);
[~, H2] = beam_shape_matrix(N2, cfg_ref.geom.L, 201);
W1 = U1*H1.';  W2 = U2*H2.';
wmax = max([max(abs(W1(:))) max(abs(W2(:)))]);
if ischar(o.exaggeration)
    f = max(0.2*L/max(wmax,eps),1);  p10 = 10^floor(log10(f));  m_ = f/p10;  ex = p10*max([1 2 5](([1 2 5]) <= m_));
else
    ex = o.exaggeration;
end
ylim_shape = 1.15*max(ex*wmax, 0.02*L);
c1 = [0 0.2 0.8];  c2 = [0.85 0.1 0.1];

% ---- figure -------------------------------------------------------------------------------------------
if inparent
    fig = o.parent;
else
    if o.show, fig = figure('name','SMA vs linear beam','numbertitle','off'); else, fig = figure('visible','off'); end
    set(fig,'units','pixels','position',[100 100 o.size],'color','w','paperpositionmode','auto');
end
ax1 = axes('parent',fig,'position',[0.12 0.56 0.84 0.36]); hold(ax1,'on')
plot(ax1,[0 L],[0 0],'k--','linewidth',1);
plot(ax1,[0 L],[0 0],'k^','markersize',7,'markerfacecolor','k');
hs2 = plot(ax1,xg, ex*W2(1,:), '-', 'color',c2, 'linewidth',2.5);
hs1 = plot(ax1,xg, ex*W1(1,:), '-', 'color',c1, 'linewidth',2.5);
axis(ax1,[0 L -ylim_shape ylim_shape]); grid(ax1,'on'); box(ax1,'on')
ylabel(ax1,sprintf('deflection x%g [mm]',ex)); xlabel(ax1,'length x [mm]');
ht = title(ax1,'');
legend(ax1,[hs1 hs2],{'SMA (loss-factor damping)','linear reference (Rayleigh)'},'location','northeast','fontsize',8);
ax2 = axes('parent',fig,'position',[0.12 0.10 0.84 0.32]); hold(ax2,'on')
sel1 = T1>=t0 & T1<=t1;  sel2 = T2>=t0 & T2<=t1;
plot(ax2,T2(sel2), ref.midspan(i2(sel2)), '-', 'color',[1 0.7 0.7]);
plot(ax2,T1(sel1), res.midspan(i1(sel1)), '-', 'color',[0.7 0.75 1]);
hm2 = plot(ax2,tf(1), m2(1), 'o', 'color',c2, 'markerfacecolor',c2, 'markersize',6);
hm1 = plot(ax2,tf(1), m1(1), 'o', 'color',c1, 'markerfacecolor',c1, 'markersize',6);
ymx = 1.1*max([res.midspan(i1(sel1)); ref.midspan(i2(sel2))]);  ymn = 1.1*min([res.midspan(i1(sel1)); ref.midspan(i2(sel2))]);
axis(ax2,[t0 t1 ymn ymx]); grid(ax2,'on'); box(ax2,'on')
xlabel(ax2,'time t [s]'); ylabel(ax2,sprintf('w(x=%g mm) [mm]', L/2));
upd = @(j) local_update(j, hs1, hs2, hm1, hm2, ex, W1, W2, m1, m2, tf, ht, speed, true);
if inparent
    upd(1);
    pl = @(j,spd) local_update(j, hs1, hs2, hm1, hm2, ex, W1, W2, m1, m2, tf, ht, spd, false);
    out = struct('gif','','nframes',nfr,'t',tf,'exaggeration',ex,'speed',speed,'fps',o.fps,'frame_size',[],'runtime_s',toc(t_start), ...
        'player',struct('name','SMA vs linear animation','tf',tf,'update',pl,'nframes',nfr));
    return
end

% ---- render ---------------------------------------------------------------------------------------------
dogif = ~isempty(o.gif);
frames = cell(nfr,1);
if dogif
    sidx = unique(round(linspace(1,nfr,min(nfr,7))));  S = [];
    for j = sidx, upd(j); S = [S; reshape(gif_grab_frame(fig), [], 3)]; end %#ok<AGROW>
    [map,lut] = gif_palette(S);
end
for j = 1:nfr
    upd(j);
    if dogif
        frames{j} = gif_quantize(gif_grab_frame(fig), lut);
    elseif o.show
        pause(0.001);
    end
    if o.verbose && mod(j,50)==0, fprintf('  comparison frame %d / %d\n', j, nfr); end
end
if ~isempty(o.still)
    upd(max(2,round(nfr/8)));  imwrite(gif_grab_frame(fig), o.still);
end
if o.show, drawnow; else, close(fig); end
out = struct('gif','','nframes',nfr,'t',tf,'exaggeration',ex,'speed',speed,'fps',o.fps,'frame_size',[],'runtime_s',0);
if dogif
    [h,w] = size(frames{1});
    gif_write(o.gif, frames, map, o.fps);
    out.gif = o.gif;  out.frame_size = [w h];
    if o.verbose
        fi = dir(o.gif);
        fprintf('comparison animation: %d frames, %gx slow motion, deflection x%g, %.1f MB -> %s\n', nfr, 1/speed, ex, fi.bytes/1e6, o.gif);
    end
end
out.runtime_s = toc(t_start);
end

function local_update(j, hs1, hs2, hm1, hm2, ex, W1, W2, m1, m2, tf, ht, speed, dodraw)
set(hs1,'ydata',ex*W1(j,:));  set(hs2,'ydata',ex*W2(j,:));
set(hm1,'xdata',tf(j),'ydata',m1(j));  set(hm2,'xdata',tf(j),'ydata',m2(j));
set(ht,'string',sprintf('t = %.4f s   (slow motion x%.3g, deflection x%g)', tf(j), 1/speed, ex),'fontsize',9);
if dodraw, drawnow; end
end
