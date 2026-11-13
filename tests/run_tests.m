function ok = run_tests(mode)
% RUN_TESTS  checks that this code still gives the same numbers as the original.
%
%   run_tests            quick tests (~3-6 min): short runs, compared with the ORIGINAL scripts
%   run_tests('full')    + the full 0.8 s base case and composite case (~10 min), peaks vs thesis
%
% What is compared: the refactored solver (run_case / sma_beam_solve) against the unmodified
% original scripts in archive/original_case_scripts/, both run with the same shortened Tf and
% the same Octave ode23t stand-in; the copy of the original is patched ONLY in the Newton loop
% condition (&& -> ||, KNOWN_ISSUES.md #1). Equal to round-off = the numerics of the main method are unchanged.
if nargin < 1, mode = 'quick'; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths(); addpath(fileparts(mfilename('fullpath')));
set(0,'defaultfigurevisible','off');
orig = fullfile('archive','original_case_scripts');
nfail = 0;

    function check(name, cond, detail)
        if nargin < 3, detail = ''; end
        if cond, fprintf('PASS  %s   %s\n', name, detail); else, fprintf('FAIL  %s   %s\n', name, detail); nfail = nfail + 1; end
    end

% 1. composite matrices vs the ones saved by hand in 2020
ref = load(fullfile(root,'tests','reference','Matrices_composite_original.mat'));
cfg = case_composite(); hle = composite_matrices(cfg.hle);
d = max([max(abs(hle.MFF_HLE(:)-ref.MFF_HLE(:))./max(abs(ref.MFF_HLE(:)))), ...
         max(abs(hle.KFF_HLE(:)-ref.KFF_HLE(:))./max(abs(ref.KFF_HLE(:)))), ...
         max(abs(hle.CFF_HLE(:)-ref.CFF_HLE(:))./max(abs(ref.CFF_HLE(:))))]);
check('composite_matrices == original Matrices.mat', d < 1e-9, sprintf('max rel diff %.2e (round-off: the 2020 file was written by MATLAB eig)', d));

% 2. base case: refactored vs original, short run
Tf = 0.03;
o = run_original('Loss_FEM1D_base.m', {'data/p.mat','p.mat'; fullfile(orig,'Loss_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf, cell(0,2));
r = run_case('base','Tf',Tf,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('base (Tf=0.03) ~= original (rel 1e-5)', same, detail);

% 2b. field figures: plot_fields writes ONE martensite figure and ONE strain/deformation figure (no empty window)
od = tempname(); mkdir(od); nfig0 = numel(findall(0,'type','figure'));
rf = run_case('base','Tf',0.03,'figures','save','plot_fields',true,'save_mat',false,'outdir',od);
check('plot_fields: <case>_fields.png, <case>_strain.png and <case>_martensite_xt.png written, no figure left open', ...
    exist(fullfile(od,'base_fields.png'),'file')==2 && exist(fullfile(od,'base_strain.png'),'file')==2 && exist(fullfile(od,'base_martensite_xt.png'),'file')==2 && numel(findall(0,'type','figure'))==nfig0, ...
    sprintf('(%d figures open)', numel(findall(0,'type','figure'))));

% 2c. mid-span value: deflection at exactly x = L/2 (Hermite) - consistent with the strain_figure shape, the
%     node value, and exact for synthetic data (src/midspan_deflection.m; dof map of build_beam.m)
Nn = 16; Ln = 1000; Lk = Ln/(Nn-1); xn = (0:Nn-1)'*Lk;
wn = sin(pi*xn/Ln); thn = pi/Ln*cos(pi*xn/Ln); un = zeros(1,2*Nn-2);
for k = 2:Nn-1, un(2*k-3) = wn(k); un(2*k-2) = thn(k); end
un(2*Nn-3) = thn(1); un(2*Nn-2) = thn(Nn);
[wmid,wnode,minfo] = midspan_deflection([zeros(1,2*Nn-2) un], Nn, Ln);
check('midspan_deflection: sin(pi x/L) -> 1 at x = L/2 (Hermite, N=16), node 9 is lower', abs(wmid-1) < 1e-4 && abs(wnode-sin(pi*xn(9)/Ln)) < 1e-12 && minfo.x_node > 500, ...
    sprintf('w(L/2) = %.7f, node %d (x = %.1f mm) = %.5f', wmid, minfo.node, minfo.x_node, wnode));
Nn = 15; Lk = Ln/(Nn-1); xn = (0:Nn-1)'*Lk; wn = sin(pi*xn/Ln); thn = pi/Ln*cos(pi*xn/Ln); un = zeros(1,2*Nn-2);
for k = 2:Nn-1, un(2*k-3) = wn(k); un(2*k-2) = thn(k); end
un(2*Nn-3) = thn(1); un(2*Nn-2) = thn(Nn);
[wmid,~,minfo] = midspan_deflection([zeros(1,2*Nn-2) un], Nn, Ln);
check('midspan_deflection: N-1 even -> the middle node itself', abs(wmid-1) < 1e-12 && minfo.x_node == 500, sprintf('node %d at x = %g mm', minfo.node, minfo.x_node));
fh = figure('visible','off');
[xs_,ws_] = strain_figure(rf.fields.epsyylon, rf.XX, 16, 1000, 10, rf.TT(end), fh); close(fh);
[~,im_] = min(abs(xs_-500));
check('res.midspan(end) == deformed shape of strain_figure at x = L/2', abs(xs_(im_)-500) < 1e-9 && abs(ws_(im_)-rf.midspan(end)) < 1e-9*max(1,abs(rf.midspan(end))), ...
    sprintf('%.6f vs %.6f mm', ws_(im_), rf.midspan(end)));
mm_ = abs(rf.midspan_node) > 0.5*max(abs(rf.midspan_node));
rat_ = rf.midspan(mm_)./rf.midspan_node(mm_);
check('mode-1 midspan is ~0.5-0.7 % larger than node 9', all(rat_ > 1.003 & rat_ < 1.010) && rf.x_midspan_node > 530 && rf.x_midspan_node < 535, ...
    sprintf('ratio %.4f..%.4f, node 9 at x = %.1f mm', min(rat_), max(rat_), rf.x_midspan_node));

% 2d. beam animation (src/animate_beam.m, beam_shape_matrix.m, strain_from_dofs.m): shape / strain consistent with the
%     model, GIF written (multi-frame, delay, looping), frame count = what the uniform frame rate says
[xg_,Hm_] = beam_shape_matrix(16, 1000, 201);
wl_ = Hm_ * rf.XX(end,2*16-1:end).';
check('beam_shape_matrix: w(L/2) == res.midspan; shape == strain_figure deformed shape', abs(wl_((numel(xg_)+1)/2) - rf.midspan(end)) < 1e-9*max(1,abs(rf.midspan(end))) && max(abs(interp1(xg_,wl_,xs_) - ws_)) < 1e-3*max(abs(ws_)), ...
    sprintf('%.6f vs %.6f mm', wl_((numel(xg_)+1)/2), rf.midspan(end)));
E_ = strain_from_dofs(rf.cfg, rf.XX(end,2*16-1:end));
check('strain_from_dofs (B*u) == res.fields.epsyylon at the last step', max(abs(E_ - rf.fields.epsyylon(:,end))) < 1e-9*max(abs(E_)), ...
    sprintf('max rel diff %.1e', max(abs(E_ - rf.fields.epsyylon(:,end)))/max(abs(E_))));
ra = run_case('base','Tf',0.07,'figures','none','save_mat',false);   % the first peak is at 0.028 s: the run is shifted to it
gf = fullfile(od,'anim_test.gif');  sf_ = fullfile(od,'anim_still.png');
ao = animate_beam(ra, ra.cfg, 'gif', gf, 'still', sf_, 'verbose', false);
inf_ = imfinfo(gf);
nexp = floor(ra.TT(end)*25/0.1) + 1;
check('animate_beam: GIF written, frame count = uniform 25 fps at 0.1x speed', exist(gf,'file')==2 && numel(inf_) == ao.nframes && ao.nframes == nexp && ao.nframes > 5, ...
    sprintf('%d frames (expected %d), %d x %d px, %.0f kB, %.1f s', numel(inf_), nexp, ao.frame_size(1), ao.frame_size(2), dir_bytes(gf)/1e3, ao.runtime_s));
check('GIF: 25 fps delay (4 cs), loops forever, LZW-compressed, still PNG written', inf_(1).DelayTime == 4 && inf_(1).LoopCount == 0 && dir_bytes(gf) < 0.15*ao.nframes*prod(ao.frame_size) && exist(sf_,'file')==2 && numel(findall(0,'type','figure'))==nfig0, ...
    sprintf('DelayTime %g cs, LoopCount %g', inf_(1).DelayTime, inf_(1).LoopCount));
rg = run_case('base','Tf',0.07,'figures','none','save_mat',false,'outdir',od,'animate','gif','animate_opts',{'speed',0.2,'strain',false});
check('run_case(''animate'',''gif''): results/<case>_animation.gif written', exist(fullfile(od,'base_animation.gif'),'file')==2 && rg.animation.nframes == floor(rg.TT(end)*25/0.2)+1 && rg.animation.nframes >= 2, ...
    sprintf('%d frames', rg.animation.nframes));
rl = run_case('hle_linear','Tf',0.07,'figures','none','save_mat',false,'outdir',od,'animate','gif','animate_opts',{'verbose',false});
check('animation of a linear case (strain from B*u, no field history)', exist(fullfile(od,'hle_linear_animation.gif'),'file')==2 && rl.animation.nframes > 5, sprintf('%d frames', rl.animation.nframes));

% 2e. looping GIFs, comparison animation, martensite snapshots / maps, dashboard
%     (a) every GIF carries the NETSCAPE loop block (gif_write); gif_ensure_loop repairs a GIF without it
fid_ = fopen(gf,'rb'); bb_ = fread(fid_,Inf,'uint8=>uint8').'; fclose(fid_);
loopblk_ = [uint8([33 255 11]) uint8('NETSCAPE2.0') uint8([3 1 0 0 0])];
check('animate_beam GIF contains the NETSCAPE loop block (loop forever)', ~isempty(strfind(char(bb_),char(loopblk_))));
gn_ = fullfile(od,'noloop.gif');  i_ = strfind(char(bb_),char(loopblk_));
fid_ = fopen(gn_,'wb'); fwrite(fid_,[bb_(1:i_(1)-1) bb_(i_(1)+numel(loopblk_):end)],'uint8'); fclose(fid_);
added_ = gif_ensure_loop(gn_);  inn_ = imfinfo(gn_);
check('gif_ensure_loop inserts the loop block into a GIF without it (still a valid GIF)', added_ && inn_(1).LoopCount == 0 && numel(inn_) == numel(inf_) && ~gif_ensure_loop(gn_));
%     (b) martensite snapshots do not change the solution; about 100 of them; shifted like TT
cs_ = case_base(); cs_.time.Tf = 0.1; cs_.post.n_snapshots = 100;
rs_ = run_case(cs_,'figures','none','save_mat',false);
r0_ = run_case('base','Tf',0.1,'figures','none','save_mat',false);
check('snapshots (cfg.post.n_snapshots) leave the solution bit-identical', isequal(rs_.XX,r0_.XX) && isequal(rs_.midspan,r0_.midspan) && ~isfield(r0_,'snapshots'));
ns_ = numel(rs_.snapshots.t);
check('snapshots: ~100 evenly spaced times within the shifted run, fractions in [0,1]', ns_ > 60 && ns_ <= 100 && all(diff(rs_.snapshots.t) > 0) && rs_.snapshots.t(1) >= 0 && rs_.snapshots.t(end) <= rs_.TT(end)+1e-12 ...
    && size(rs_.snapshots.xisp,1) == 750 && min([rs_.snapshots.xisp(:); rs_.snapshots.xism(:)]) >= 0 && max(rs_.snapshots.xiT(:)) <= 1+1e-12, sprintf('%d snapshots', ns_));
rf2_ = run_case('base','Tf',0.1,'figures','none','save_mat',false,'plot_fields',true);
check('last snapshot == fields of the last step (xi_s+, xi_s-, xi_T)', isfield(rf2_,'snapshots') && abs(rf2_.snapshots.t(end) - rf2_.TT(end)) < 0.03*rf2_.TT(end) ...
    && max(abs(rf2_.snapshots.xisp(:,end) - rf2_.fields.xysp(:,end))) < 0.05*max(rf2_.fields.xysp(:,end)));
%     (c) static x-t map and the fields animation
fx_ = figure('visible','off'); martensite_xt(rs_, cs_, fx_); fnx_ = fullfile(od,'xt_test.png'); print(fx_,'-dpng','-r60',fnx_); close(fx_);
Mxt_ = snapshot_fields(rs_,cs_,'xisp','mean');
check('martensite_xt: figure drawn and saved; section average of xi_s+ is non-negative, grows in time at mid-span', exist(fnx_,'file')==2 && all(Mxt_(:) >= 0) && Mxt_(end,round(end/2)) > Mxt_(1,round(end/2)), sprintf('mid-span %.2e -> %.2e', Mxt_(1,round(end/2)), Mxt_(end,round(end/2))));
gfl_ = fullfile(od,'fields_test.gif');
af_ = animate_fields(rf2_, rf2_.cfg, 'gif', gfl_, 'step', 3, 'verbose', false);
inf2_ = imfinfo(gfl_);
check('animate_fields: GIF written with one frame per (3rd) snapshot, loops forever', numel(inf2_) == af_.nframes && af_.nframes == numel(1:3:numel(rf2_.snapshots.t)) && inf2_(1).LoopCount == 0, sprintf('%d frames, %.0f kB', af_.nframes, dir_bytes(gfl_)/1e3));
%     (d) SMA vs linear comparison animation
rc_ = run_case('base','Tf',0.07,'figures','none','save_mat',false,'outdir',od,'animate_compare','gif','animate_opts',{'speed',0.2,'verbose',false});
inf3_ = imfinfo(fullfile(od,'base_animation_compare.gif'));
check('run_case(''animate_compare'',''gif''): results/<case>_animation_compare.gif written, loops, linear reference run', numel(inf3_) == rc_.animation_compare.nframes && inf3_(1).LoopCount == 0 && isfield(rc_,'reference') && rc_.animation_compare.nframes > 5, sprintf('%d frames', numel(inf3_)));
%     (e) dashboard: tab buttons and panels built, switching works (needs uipanel: Octave Qt / MATLAB)
if dashboard_supported()
    fdb_ = dashboard(rf2_, rf2_.cfg, [], 'visible', false);
    dd_ = getappdata(fdb_,'dash');
    dashboard_select(fdb_,numel(dd_.names));
    check('dashboard: tabs Deflection / Strain / Martensite, only the selected panel visible', numel(dd_.names) == 4 && strcmp(get(dd_.panels(end),'visible'),'on') && strcmp(get(dd_.panels(1),'visible'),'off') && get(dd_.buttons(end),'value') == 1, strjoin(dd_.names,' | '));
    close(fdb_);
else
    fprintf('SKIP  dashboard (graphics toolkit without uipanel)\n');
end

% 2f. one-window animations: tabs, continuous loop, speed control, pause / stop / closing the window
if dashboard_supported()
    rp_ = run_case('base','Tf',0.1,'figures','none','save_mat',false,'plot_fields',true,'animate_compare','gif','outdir',od,'animate_opts',{'speed',0.5,'verbose',false});
    A_ = struct('beam',true,'compare',true,'fields',true,'refc',rp_.reference,'speed',0.5,'autoplay',true);
    nfig0_ = numel(findobj(0,'type','figure'));
    fp_ = dashboard(rp_, rp_.cfg, [], 'visible', false, 'plots', false, 'anim', A_);
    dp_ = getappdata(fp_,'dash');
    check('dashboard: animation tabs Beam / SMA vs linear / Martensite, no other window opened', isequal(dp_.names,{'Beam animation','SMA vs linear animation','Martensite animation'}) && all(dp_.isanim) && numel(findobj(0,'type','figure')) == nfig0_+1);
    ii_ = dashboard_play(fp_,'maxtime',1.5);
    check('dashboard_play: only the visible tab is animated, the hidden ones wait', ii_.frames(1) > 5 && all(ii_.frames(2:3) == 0) && ~ii_.closed, sprintf('frames per tab %s in 1.5 s', mat2str(ii_.frames)));
    dashboard_select(fp_,3);  ii_ = dashboard_play(fp_,'maxtime',1.0);
    check('switching the tab moves the animation to the new tab', ii_.frames(3) > 3 && all(ii_.frames(1:2) == 0), sprintf('frames per tab %s', mat2str(ii_.frames)));
    %   continuous loop: wraps around at the end; speed control scales the playback
    dashboard_select(fp_,1);
    dashboard_ctrl(fp_,1,'speed',1);   i1_ = dashboard_play(fp_,'maxtime',1.5);
    dashboard_ctrl(fp_,1,'speed',0.25);  i2_ = dashboard_play(fp_,'maxtime',1.5);
    check('animation loops continuously and the speed control scales it (1x vs 0.25x)', i1_.wraps(1) >= 5 && i2_.wraps(1) >= 1 && i1_.wraps(1) > 2*i2_.wraps(1), sprintf('passes in 1.5 s: %d at 1x, %d at 0.25x (Tf = %g s)', i1_.wraps(1), i2_.wraps(1), rp_.TT(end)));
    dashboard_ctrl(fp_,1,'slower'); sp_ = getappdata(fp_,'play_speed'); dashboard_ctrl(fp_,1,'faster'); sp2_ = getappdata(fp_,'play_speed');
    for q_ = 1:20, dashboard_ctrl(fp_,1,'slower'); end;  spm_ = getappdata(fp_,'play_speed');
    for q_ = 1:30, dashboard_ctrl(fp_,1,'faster'); end;  spM_ = getappdata(fp_,'play_speed');
    hc_ = getappdata(fp_,'play_ctrl');
    check('speed buttons: -/+ step by 1.5x and stay within 0.02 ... 1; the slider and the label follow', abs(sp_(1)-0.25/1.5) < 1e-12 && abs(sp2_(1)-0.25) < 1e-12 && abs(spm_(1)-0.02) < 1e-12 && abs(spM_(1)-1) < 1e-12 && abs(get(hc_(1,2),'value')) < 1e-12 && ~isempty(strfind(get(hc_(1,3),'string'),'1x real time')), sprintf('%.4g %.4g %.4g %.4g', sp_(1), sp2_(1), spm_(1), spM_(1)));
    set(hc_(1,2),'value',log10(0.1)); dashboard_ctrl(fp_,1,'slider'); sp3_ = getappdata(fp_,'play_speed');
    check('speed slider sets the speed (log scale)', abs(sp3_(1)-0.1) < 1e-9);
    %   pause / Stop / closing the window while the loop runs
    dashboard_ctrl(fp_,1,'speed',1);
    ip_ = dashboard_play(fp_,'maxtime',1.5,'hook',@(f,i) play_hook(f,i,'pause',5));
    onp_ = getappdata(fp_,'play_on');
    check('Pause: the tab stops advancing (loop keeps running, button reads Play)', ip_.frames(1) > 0 && ip_.frames(1) <= 8 && ~onp_(1) && strcmp(get(hc_(1,1),'string'),'Play') && ip_.elapsed >= 1.4, sprintf('%d frames before the pause, play_on %d, button %s, %.2f s', ip_.frames(1), onp_(1), get(hc_(1,1),'string'), ip_.elapsed));
    close(fp_);
    fp_ = dashboard(rp_, rp_.cfg, [], 'visible', false, 'plots', false, 'anim', A_);
    is_ = dashboard_play(fp_,'maxtime',5,'hook',@(f,i) play_hook(f,i,'stop',10));
    check('Stop ends the loop at once and leaves the window open', is_.elapsed < 2 && ishandle(fp_) && ~getappdata(fp_,'play_run') && ~any(getappdata(fp_,'play_on')), sprintf('left after %.2f s', is_.elapsed));
    setappdata(fp_,'play_on',true(1,3));  ir_ = dashboard_play(fp_,'maxtime',1);
    check('the loop can be started again after Stop (Play)', ir_.frames(1) > 5, sprintf('%d frames', ir_.frames(1)));
    fp2_ = dashboard(rp_, rp_.cfg, [], 'visible', false, 'plots', false, 'anim', A_);
    ok_ = true;  try, txt_ = evalc('ic_ = dashboard_play(fp2_,''maxtime'',5,''hook'',@(f,i) play_hook(f,i,''close'',10));'); catch, ok_ = false; end
    check('closing the window during playback ends the loop without an error', ok_ && ic_.closed && ic_.elapsed < 2 && ~ishandle(fp2_) && isempty(strtrim(txt_)), sprintf('left after %.2f s', ic_.elapsed));
    if ishandle(fp_), close(fp_); end
    %   run_case: animations in the dashboard, GIF speed, blocking behaviour
    nf0_ = numel(findobj(0,'type','figure'));
    t0_ = tic;
    setappdata(0,'sma_play_test',struct('maxtime',1));
    rq_ = run_case('base','Tf',0.05,'figures','save','save_mat',false,'outdir',od,'animate',true,'animate_compare','both','animate_fields',true,'plot_fields',true,'anim_speed',0.5);
    tq_ = toc(t0_);  rmappdata(0,'sma_play_test');
    dq_ = getappdata(rq_.dashboard,'dash');
    check('run_case(animate*): animations are tabs of ONE window, no extra windows, loops until stopped (here: 1 s)', isfield(rq_,'dashboard') && ishandle(rq_.dashboard) && numel(findobj(0,'type','figure')) == nf0_+1 && numel(dq_.names) == 3 && ~isfield(rq_,'animation') && isfield(rq_,'animation_compare'), strjoin(dq_.names,' | '));
    close(rq_.dashboard);
    setappdata(0,'sma_play_test',struct('hook',@(f,i) play_hook(f,i,'close',15)));
    t0_ = tic;  rq2_ = run_case('base','Tf',0.05,'figures','none','save_mat',false,'animate',true);  tq2_ = toc(t0_);  rmappdata(0,'sma_play_test');
    check('run_case returns once the animation window is closed', tq2_ < 15 && ~ishandle(rq2_.dashboard), sprintf('%.1f s in all (solve + build + playback)', tq2_));
    rq3_ = run_case('base','Tf',0.05,'figures','none','save_mat',false,'animate',true,'anim_block',false);
    onq_ = getappdata(rq3_.dashboard,'play_on');
    check('anim_block = false: returns at once with the animations paused (Play starts them)', ishandle(rq3_.dashboard) && ~any(onq_));
    close(rq3_.dashboard);
    ga_ = run_case('base','Tf',0.07,'figures','none','save_mat',false,'outdir',od,'animate_compare','gif','anim_speed',0.5);  ia_ = imfinfo(fullfile(od,'base_animation_compare.gif'));
    gb_ = run_case('base','Tf',0.07,'figures','none','save_mat',false,'outdir',od,'animate_compare','gif','anim_speed',0.25);  ib_ = imfinfo(fullfile(od,'base_animation_compare.gif'));
    check('anim_speed sets the GIF speed (half the speed = twice the frames)', abs(numel(ib_) - 2*numel(ia_)) <= 2 && abs(ga_.animation_compare.speed-0.5) < 1e-12, sprintf('%d frames at 0.5, %d at 0.25', numel(ia_), numel(ib_)));
    gf1_ = run_case('base','Tf',0.1,'figures','none','save_mat',false,'outdir',od,'animate_fields','gif','anim_speed',0.2,'plot_fields',true);
    gf2_ = run_case('base','Tf',0.1,'figures','none','save_mat',false,'outdir',od,'animate_fields','gif','anim_speed',0.1,'plot_fields',true);
    check('anim_speed for the martensite GIF: frames at uniform times, scale with 1/speed', gf2_.animation_fields.nframes > 1.7*gf1_.animation_fields.nframes, sprintf('%d vs %d frames', gf1_.animation_fields.nframes, gf2_.animation_fields.nframes));
    ok_ = false; try, run_case('base','Tf',0.01,'figures','none','save_mat',false,'anim_speed',-1); catch, ok_ = true; end
    check('anim_speed must be a positive number', ok_);
else
    fprintf('SKIP  one-window animations (graphics toolkit without uipanel)\n');
end
%     colour bar tick labels: 5 evenly spaced, increasing ticks from clim(1) to clim(2), each label is its tick value
fcb_ = figure('visible','off');
okcb_ = true;  dcb_ = '';
for cl_ = {[0 0.0126], [0.9744 1], [-0.31e-3 0.31e-3], [0 1e-12]}
    axb_ = hand_colorbar(fcb_,[0.8 0.1 0.03 0.8],cl_{1},jet(128),'','%.4f');
    tk_ = get(axb_,'ytick');  lb_ = str2double(get(axb_,'yticklabel'));  yl_ = get(axb_,'ylim');
    good_ = numel(tk_) == 5 && all(diff(tk_) > 0) && max(abs(diff(tk_)-mean(diff(tk_)))) < 1e-9*max(1,abs(cl_{1}(2))) ...
        && abs(tk_(1)-cl_{1}(1)) < 1e-12 && abs(tk_(end)-cl_{1}(2)) < 1e-12 && all(isfinite(lb_)) && all(diff(lb_(:).') >= 0) ...
        && all(abs(lb_(:).'-tk_) <= 0.5e-4 + 1e-12) && tk_(1) >= yl_(1)-1e-12 && tk_(end) <= yl_(2)+1e-12;
    okcb_ = okcb_ && good_;  dcb_ = [dcb_ sprintf('[%g %g]: %s  ', cl_{1}, strjoin(get(axb_,'yticklabel')',' '))];
end
check('hand_colorbar: 5 evenly spaced, increasing tick labels equal to the tick values, ends at the colour limits', okcb_, dcb_);
close(fcb_);
if dashboard_supported()
    fcf_ = figure('visible','off');  afo_ = animate_fields(rp_, rp_.cfg, 'parent', fcf_, 'verbose', false);
    cbs_ = findobj(fcf_,'type','axes','-and','yaxislocation','right');  okf_ = numel(cbs_) == 3;  dfl_ = '';
    for q_ = 1:numel(cbs_)
        tk_ = get(cbs_(q_),'ytick');  lb_ = str2double(get(cbs_(q_),'yticklabel'));
        okf_ = okf_ && all(diff(tk_) > 0) && all(diff(lb_(:).') >= 0) && all(abs(lb_(:).'-tk_) <= 0.5e-4+1e-12) && any(all(abs(afo_.clim - repmat([tk_(1) tk_(end)],3,1)) < 1e-12, 2));
        dfl_ = [dfl_ strjoin(get(cbs_(q_),'yticklabel')',' ') ' | '];
    end
    check('animate_fields colour bars: tick labels increasing, equal to the ticks, ends = colour limits of each map', okf_, dfl_);
    close(fcf_);
    fcb2_ = figure('visible','off');  abo_ = animate_beam(rf, rf.cfg, 'parent', fcb2_, 'verbose', false);
    cb2_ = findobj(fcb2_,'type','axes','-and','yaxislocation','right');
    tk_ = get(cb2_,'ytick');  lb_ = str2double(get(cb2_,'yticklabel'));
    check('animate_beam strain colour bar: labels = tick values x 1e3, increasing', numel(cb2_) == 1 && all(diff(tk_) > 0) && all(diff(lb_(:).') > 0) && all(abs(lb_(:).'-tk_*1e3) <= 0.5e-2+1e-12), strjoin(get(cb2_,'yticklabel')',' '));
    close(fcb2_);
end
%     composite: the martensite / strain fields exist only over the SMA core (geom.thick), the GFRP skins are linear
cb_ = case_composite();  cs0_ = case_base();
check('section label: SMA core named for the sandwich beam, plain h for the SMA beam', ~isempty(strfind(sma_section_label(cb_,'axis'),'SMA core')) && ~isempty(strfind(sma_section_label(cb_,'suffix'),'pm5')) && strcmp(sma_section_label(cs0_,'axis'),'h [mm]') && isempty(sma_section_label(cs0_,'suffix')) && cb_.geom.thick == 10);
rcf_ = run_case('composite','Tf',0.02,'figures','none','save_mat',false,'plot_fields',true);
check('composite: martensite snapshots are defined at the SMA integration points only (n_y x n_csi x elements)', size(rcf_.snapshots.xisp,1) == rcf_.snapshots.n_y*rcf_.snapshots.n_csi*(cb_.geom.N-1) && cb_.hle.thick > cb_.hle.thick2 && cb_.hle.thick2 == cb_.geom.thick, sprintf('%d points, SMA core %g mm of %g mm', size(rcf_.snapshots.xisp,1), cb_.geom.thick, cb_.hle.thick));
fm_ = figure('visible','off');  ff_ = rcf_.fields;
meshare(ff_.epsyylon,ff_.xysp,ff_.xysm,ff_.xyT,rcf_.TT,cb_.geom.thick,cb_.geom.L,fm_,struct('ylabel',sma_section_label(cb_,'axis'),'suffix',sma_section_label(cb_,'suffix')));
axm_ = findobj(fm_,'type','axes');  axm_ = axm_(arrayfun(@(a) ~strcmp(get(a,'tag'),'colorbar') && ~isempty(get(get(a,'ylabel'),'string')) && ~isempty(strfind(get(get(a,'ylabel'),'string'),'h [mm]')), axm_));
yl_ = cell2mat(arrayfun(@(a) get(a,'ylim'), axm_(:), 'uniformoutput', false));
check('martensite field plot: h axis spans exactly the SMA core +-thick/2, labelled SMA core', ~isempty(axm_) && all(abs(yl_(:,1)+cb_.geom.thick/2) < 1e-12) && all(abs(yl_(:,2)-cb_.geom.thick/2) < 1e-12), sprintf('ylim %s', mat2str(yl_(1,:))));
close(fm_);


% 3. global state is reset between runs: same case twice in one session -> identical
r2 = run_case('base','Tf',Tf,'figures','none','save_mat',false);
check('re-run in the same session identical', isequal(r.XX,r2.XX), sprintf('max |dXX| = %g', max(abs(r.XX(:)-r2.XX(:)))));

% 4. composite
o = run_original('Loss_Composite_FEM1D.m', {'data/p.mat','p.mat'; 'tests/reference/Matrices_composite_original.mat','Matrices.mat'; fullfile(orig,'Loss_Composite_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf, cell(0,2));
% (1) with the 2020 Matrices.mat the result is bit-identical to the original script
cfgc = case_composite(); cfgc.time.Tf = Tf;
cfgc.hle.matrices_file = fullfile(root,'tests','reference','Matrices_composite_original.mat');
r = run_case(cfgc,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('composite (Tf=0.03, original Matrices.mat) ~= original (rel 1e-5)', same, detail);
% (2) with the matrices recomputed by composite_matrices.m: ~1e-12 round-off in the Rayleigh coefficients
%     (the 2020 file was written by MATLAB), which gives ~1e-9 relative differences in the result
r = run_case('composite','Tf',Tf,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('composite (Tf=0.03, recomputed matrices) ~= original', same, detail);

% 5. forced pulse (Loss/Transient): short run, includes the ramp and the 0.02 s release
Tf5 = 0.025;
o = run_original('Loss_Transient_FEM1D.m', {'data/p.mat','p.mat'; fullfile(orig,'Loss_Transient_myfun.m'),'myfun.m'; fullfile(orig,'Loss_Transient_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf5, cell(0,2));
r = run_case('pulse','Tf',Tf5,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('pulse (Tf=0.025) ~= original (rel 1e-5)', same, detail);

% 6. super-elastic pre-strain (ode15s of the original replaced by the ode23t stand-in, as in the port)
Tf6 = 0.02;
o = run_original('Loss_SE_pre-strain_FEM1D.m', {'data/p_SE.mat','p.mat'; fullfile(orig,'Loss_SE_pre-strain_myfun.m'),'myfun.m'; fullfile(orig,'Loss_fun_ode_1D.m'),'fun_ode_1D.m'}, Tf6, {'ode15s(@','ode23t(@'});
r = run_case('prestrain_se','Tf',Tf6,'figures','none','save_mat',false);
[same,detail] = compare_runs(o,r); check('prestrain_se (Tf=0.02) ~= original (rel 1e-5)', same, detail);

if strcmp(mode,'full')
    % 7. full base case: peaks vs thesis Fig. 4.3 and the stored reference run of the Octave port
    r = run_case('base','figures','save','save_mat',false);
    pk = r.peaks(:,2)';
    expect = [7.13 5.52 4.17 3.09 2.25 1.62];   % deflection at x = L/2 (Hermite), corrected (||) Newton loop; node 9 (x = 533 mm, used by earlier versions) gives 7.09 5.49 4.15 3.07 2.24 1.61
    check('base full run: ring-down peaks', numel(pk) >= 6 && all(abs(pk(1:6)-expect) < 0.02), ...
        sprintf('peaks [mm]: %s (expected %s; thesis Fig. 4.3: ~5.5 4.15 3.1 2.2 1.6 after the first)', mat2str(pk(1:min(6,end)),3), mat2str(expect)));
    r = run_case('composite','figures','save','save_mat',false);
    check('composite full run finished', r.TT(end) > 0.7, sprintf('t_end = %.3f s, %d steps, %.0f s', r.TT(end), numel(r.TT), r.runtime_s));
end
ok = (nfail == 0);
if ok, fprintf('\nALL TESTS PASSED\n'); else, fprintf('\n%d TEST(S) FAILED\n', nfail); end
end

function b = dir_bytes(f)
d = dir(f); b = d.bytes;
end

function [same,detail] = compare_runs(o,r,tol)
if nargin < 3, tol = 1e-5; end   % relative to max|X|; round-off level (the fast path sums in a different order, see PERF_NOTES.md)
% o: original (TT,XX); r: refactored result
if numel(o.TT) ~= numel(r.TT)
    same = false; detail = sprintf('different number of steps: original %d, new %d', numel(o.TT), numel(r.TT)); return
end
dX = max(abs(o.XX(:)-r.XX(:))); dT = max(abs(o.TT(:)-r.TT(:)));
scale = max(abs(o.XX(:)));
same = dX <= tol*scale && dT <= tol*max(abs(o.TT(:)));
detail = sprintf('%d steps, max|dXX| = %.2e (rel %.1e), max|dT| = %.1e', numel(r.TT), dX, dX/scale, dT);
end
