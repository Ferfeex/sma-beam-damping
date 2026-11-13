function out = animate_beam(res, cfg, varargin)
% ANIMATE_BEAM  post-run animation of the vibrating beam (and optional GIF export).
%
%   out = animate_beam(res, cfg, 'name', value, ...)       res = result of run_case / sma_beam_solve
%
%   Panels (top to bottom)
%     1. deflected beam, cubic Hermite shape (nodal w and theta, src/beam_shape_matrix.m), against the
%        undeformed axis, deflection exaggerated by a factor shown in the label
%     2. (option 'strain') strain field over length x thickness, fixed colour limits over the whole run
%        (strain recomputed from the nodal dofs with the B matrix, no field history needed)
%     3. (option 'history') mid-span deflection (x = L/2) versus time with a moving marker
%
%   Options (defaults)
%     'gif'          ''     file name of the GIF to write ('' = no file)
%     'show'         false  play in a visible window (the frames are shown while they are rendered)
%     'fps'          25     frames per second of the animation / GIF
%     'speed'        'auto' playback speed, simulated seconds per real second. 'auto' = 0.1 (10x slow motion)
%                           for an oscillation of ~10 Hz; a number gives 1 = real time, 0.1 = ten times slower
%     'exaggeration' 'auto' multiplication of the deflection in panel 1 (number >= 1)
%     'strain'       true   strain panel
%     'history'      true   mid-span time-history panel with moving marker
%     'tspan'        []     [t0 t1] part of the run to animate (default: whole run)
%     'maxframes'    600    safety cap on the number of frames
%     'size'         [560 480]  frame size in pixels
%     'still'        ''     file name of a PNG with one still frame (the frame of largest downward mid-span
%                           deflection in the first fifth of the run)
%     'verbose'      true
%     'parent'       []     figure or uipanel to build the animation in (the dashboard tab) instead of a window or a GIF:
%                           nothing is played or written; out.player holds the frame update function (see dashboard_play.m)
%     'playframes'   300    number of uniform frame times prepared for an on-screen player ('parent' mode)
%   Returns out: .gif (file), .nframes, .t (frame times), .exaggeration, .speed, .fps, .runtime_s, .frame_size
%
%   GIF writing is pure Octave/MATLAB (src/gif_write.m: imwrite 'gif', 'DelayTime', 'LoopCount' = Inf, so the
%   GIF loops forever); the frames are quantised
%   to one common 256-colour palette. No external tool is needed. In Octave imwrite uses GraphicsMagick
%   (shipped with the Windows installer).
t_start = tic;
o = struct('gif','','show',false,'fps',25,'speed','auto','exaggeration','auto','strain',true,'history',true, ...
           'tspan',[],'maxframes',600,'size',[560 480],'still','','verbose',true,'parent',[],'playframes',300);
for k = 1:2:numel(varargin)
    if ~isfield(o,varargin{k}), error('animate_beam: unknown option "%s"', varargin{k}); end
    o.(varargin{k}) = varargin{k+1};
end
N = cfg.geom.N;  L = cfg.geom.L;  thick = cfg.geom.thick;
TT = res.TT(:);  XX = res.XX;
[TT,iu] = unique(TT);  XX = XX(iu,:);                % strictly increasing time axis for interp1
if isempty(o.tspan), o.tspan = [TT(1) TT(end)]; end
t0 = max(o.tspan(1),TT(1));  t1 = min(o.tspan(2),TT(end));
if t1 <= t0, error('animate_beam: empty time span'); end

% ---- uniform frame times ----------------------------------------------------------------------
inparent = ~isempty(o.parent);
if ischar(o.speed), speed = 0.1; else, speed = o.speed; end
nfr = max(2, floor((t1-t0)*o.fps/speed) + 1);
if inparent, nfr = o.playframes; end                  % player: frames are picked from the wall clock, not from a frame rate
if ~inparent && nfr > o.maxframes
    if o.verbose, fprintf('animate_beam: %d frames needed, capped at %d (playback speed raised to %.3g)\n', nfr, o.maxframes, (t1-t0)*o.fps/(o.maxframes-1)); end
    nfr = o.maxframes;  speed = (t1-t0)*o.fps/(nfr-1);
end
tf = linspace(t0, t1, nfr).';
U  = interp1(TT, XX(:,2*N-1:end), tf, 'linear');      % nfr x (2N-2) nodal displacements at the frame times
mid = interp1(TT, res.midspan(iu), tf, 'linear');

% ---- shape ---------------------------------------------------------------------------------------
[xg, Hm] = beam_shape_matrix(N, L, 201);
W = U * Hm.';                                          % nfr x nx deflection [mm]
wmax = max(abs(W(:)));
if ischar(o.exaggeration)
    f = 0.2*L/max(wmax,eps);                           % peak deflection drawn as ~20 % of the length
    f = max(f,1);  p10 = 10^floor(log10(f));  m_ = f/p10;  ex = p10*max([1 2 5](([1 2 5]) <= m_));   % 1, 2 or 5 x 10^k
else
    ex = o.exaggeration;
end
ylim_shape = 1.15*max(ex*wmax, 0.02*L);

% ---- strain field (fixed colour limits) ---------------------------------------------------------
if o.strain
    EPS = strain_from_dofs(cfg, U);                    % nip x nfr, like res.fields.epsyylon
    ny = cfg.geom.n_y;  nc = size(EPS,1)/ny;
    xs_ = linspace(0,L,nc).';  ys_ = linspace(thick/2,-thick/2,ny).';
    [Xs,Ys] = meshgrid(xs_,ys_);
    [xi,yi] = meshgrid(linspace(0,L,240), linspace(thick/2,-thick/2,60));
    smoothfield = @(e) interp2(Xs,Ys,smooth_mean2d(reshape(e,ny,nc),1,8),xi,yi);
    Emin = Inf; Emax = -Inf;
    for j = 1:nfr
        sm = smooth_mean2d(reshape(EPS(:,j),ny,nc),1,8);
        Emin = min(Emin,min(sm(:))); Emax = max(Emax,max(sm(:)));
    end
    if Emin < 0 && Emax > 0, Ea_ = max(abs([Emin Emax])); clim = [-Ea_ Ea_]; else, clim = [Emin Emax]; end   % symmetric when the strain changes sign (SE pre-strain: stays positive)
    if diff(clim) <= 0, clim = clim + [-1 1]*1e-12; end
end

% ---- figure ---------------------------------------------------------------------------------------
if inparent
    fig = o.parent;  par = fig;
else
    if o.show, fig = figure('name','beam animation','numbertitle','off'); else, fig = figure('visible','off'); end
    set(fig,'units','pixels','position',[100 100 o.size],'color','w','paperpositionmode','auto');
    par = fig;
end
np = 1 + o.strain + o.history;
gap = 0.055; top = 0.93; bot = 0.10; h_each = (top-bot-gap*(np-1))/np;
ypos = @(i) top - i*h_each - (i-1)*gap;             % bottom edge of panel i (1-based)
% panel 1
ax1 = axes('parent',par,'position',[0.12 ypos(1) 0.74 h_each]); hold(ax1,'on')
plot(ax1,[0 L],[0 0],'k--','linewidth',1);
plot(ax1,[0 L],[0 0],'k^','markersize',7,'markerfacecolor','k');
hs = plot(ax1,xg, ex*W(1,:), 'b', 'linewidth', 2.5);
axis(ax1,[0 L -ylim_shape ylim_shape]); grid(ax1,'on'); box(ax1,'on')
ylabel(ax1,sprintf('deflection x%g [mm]',ex)); 
set(ax1,'xticklabel',[]);
ht = title(ax1,'');
% panel 2
if o.strain
    ax2 = axes('parent',par,'position',[0.12 ypos(2) 0.74 h_each]);
    him = imagesc(ax2,xi(1,:), yi(:,1), smoothfield(EPS(:,1))); set(ax2,'ydir','normal'); set(ax2,'clim',clim); colormap(ax2,jet(128));
    ylabel(ax2,sma_section_label(cfg,'thickness')); title(ax2,['strain field' sma_section_label(cfg,'suffix')],'fontsize',9);
    if o.history, set(ax2,'xticklabel',[]); else, xlabel(ax2,'length [mm]'); end
    % colour bar drawn by hand (fixed position, same in Octave and MATLAB)
    axc = axes('parent',par,'position',[0.88 ypos(2) 0.025 h_each]);
    imagesc(axc,[0 1], linspace(clim(1),clim(2),128), linspace(clim(1),clim(2),128).'); set(axc,'ydir','normal');
    set(axc,'clim',clim); colormap(axc,jet(128)); tk = linspace(clim(1),clim(2),5);
    set(axc,'xtick',[],'yaxislocation','right','fontsize',9,'ytick',tk,'yticklabel',arrayfun(@(v) sprintf('%.2f',v*1e3),tk,'uniformoutput',false));
    title(axc,'\epsilon x 10^{-3}','fontsize',8);
end
% panel 3
if o.history
    ax3 = axes('parent',par,'position',[0.12 ypos(np) 0.74 h_each]); hold(ax3,'on')
    plot(ax3,TT(TT>=t0 & TT<=t1), res.midspan(iu(TT>=t0 & TT<=t1)), 'color',[0.6 0.6 0.6]);
    hm = plot(ax3,tf(1), mid(1), 'ro', 'markerfacecolor','r', 'markersize',6);
    axis(ax3,[t0 t1 1.1*min(res.midspan) 1.1*max(res.midspan)]); grid(ax3,'on'); box(ax3,'on')
    xlabel(ax3,'time t [s]'); ylabel(ax3,sprintf('w(x=%g mm) [mm]', L/2));
end

% ---- palette from a few sample frames, then render every frame ----------------------------------
if o.strain, him_h = him; sf = smoothfield; EP = EPS; else, him_h = []; sf = []; EP = []; end
if o.history, hm_h = hm; else, hm_h = []; end
upd = @(j) local_update(j, hs, ex, W, ht, tf, speed, him_h, sf, EP, hm_h, mid, true);
if inparent
    upd(1);
    pl = @(j,spd) local_update(j, hs, ex, W, ht, tf, spd, him_h, sf, EP, hm_h, mid, false);
    out = struct('gif','','nframes',nfr,'t',tf,'exaggeration',ex,'speed',speed,'fps',o.fps,'frame_size',[],'runtime_s',toc(t_start), ...
        'player',struct('name','Beam animation','tf',tf,'update',pl,'nframes',nfr));
    return
end
dogif = ~isempty(o.gif);
frames = cell(nfr,1);
if dogif
    ns = min(nfr,7); sidx = unique(round(linspace(1,nfr,ns)));
    S = [];
    for j = sidx
        upd(j); S = [S; reshape(gif_grab_frame(fig), [], 3)]; %#ok<AGROW>
    end
    [map, lut] = gif_palette(S);
end
tprint = tic;
for j = 1:nfr
    upd(j);
    if dogif
        A = gif_grab_frame(fig);
        key = bitshift(uint32(A(:,:,1)),-3)*1024 + bitshift(uint32(A(:,:,2)),-3)*32 + bitshift(uint32(A(:,:,3)),-3);
        frames{j} = reshape(lut(double(key)+1), size(A,1), size(A,2));
    elseif o.show
        pause(0.001);
    end
    if o.verbose && mod(j,50)==0, fprintf('  animation frame %d / %d\n', j, nfr); end
end
if ~isempty(o.still)
    [~,js] = min(mid(1:max(2,ceil(nfr/5))));  % nfr >= 2
    upd(js);  imwrite(gif_grab_frame(fig), o.still);
end
if o.show, drawnow; else, close(fig); end

out = struct('gif','','nframes',nfr,'t',tf,'exaggeration',ex,'speed',speed,'fps',o.fps,'frame_size',[],'runtime_s',0);
if dogif
    [h,w] = size(frames{1});
    d = dirname(o.gif); if ~isempty(d) && ~exist(d,'dir'), mkdir(d); end
    gif_write(o.gif, frames, map, o.fps);              % looping GIF (see gif_write.m)
    out.gif = o.gif;  out.frame_size = [w h];
    if o.verbose
        fi = dir(o.gif);
        fprintf('animation: %d frames, %gx slow motion (playback %.3g sim-s per s), deflection x%g, %.1f MB -> %s\n', ...
            nfr, 1/speed, speed, ex, fi.bytes/1e6, o.gif);
    end
end
out.runtime_s = toc(t_start);
end

% ======================================================================================================
function local_update(j, hs, ex, W, ht, tf, speed, him, sf, EP, hm, mid, dodraw)
set(hs,'ydata',ex*W(j,:));
set(ht,'string',sprintf('t = %.4f s   (slow motion x%.3g, deflection x%g)', tf(j), 1/speed, ex),'fontsize',9);
if ~isempty(him), set(him,'cdata',sf(EP(:,j))); end
if ~isempty(hm), set(hm,'xdata',tf(j),'ydata',mid(j)); end
if dodraw, drawnow; end
end

function d = dirname(f)
[d,~,~] = fileparts(f);
end
