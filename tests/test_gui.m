function test_gui(shotdir)
% TEST_GUI  automatic test of the start window (needs a display, e.g. xvfb-run on Linux):
% opens sma_gui, sets values programmatically, triggers the Run callback with a short Tf.
%   xvfb-run -a octave --no-gui tests/test_gui.m          (or call test_gui('<folder for screenshots>'))
if nargin < 1, shotdir = ''; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
nfail = 0;
    function check(name,ok,detail)
        if nargin < 3, detail = ''; end
        if ok, fprintf('PASS  %s %s\n',name,detail); else, fprintf('FAIL  %s %s\n',name,detail); nfail = nfail+1; end
    end
% reference: scripted run without the window
ref = run_case('base','Tf',0.05,'figures','none','save_mat',false);
pfile_before = dir(fullfile(root,'data','p.mat'));
d0 = load(fullfile(root,'data','p.mat'));

setappdata(0,'sma_play_test',struct('maxtime',1));    % animations normally loop until stopped; here each playback ends after 1 s
fig = sma_gui();  setappdata(fig,'Silent',true);  drawnow;
S = getappdata(fig,'S'); H = S.H;
check('window opens', ishandle(fig) && numel(S.C)==7, sprintf('(%d cases)',numel(S.C)));
check('p fields pre-filled from data/p.mat', all(~cellfun(@isempty,get(H.pe,'string'))), ...
    sprintf('[%s]', strjoin(get(H.pe,'string').',' ')));
check('runtime estimate shown', ~isempty(get(H.est,'string')));

% --- every case can be selected; linear cases disable p / compare ----------------------------------
okall = true;
for k = 1:numel(S.C)
    set(H.list,'value',k); sma_gui('select',fig);
    lin = strcmp(S.C(k).model,'linear');
    en = get(H.pe(1),'enable');
    okall = okall && strcmp(en, ternary(lin,'off','on')) && strcmp(get(H.cmp,'enable'), ternary(lin,'off','on'));
end
check('case selection toggles p fields / compare for linear cases', okall);
set(H.list,'value',4); sma_gui('select',fig);
check('pulse case: estimate shown, no "very slow" warning with the fast code path', ~isempty(get(H.est,'string')) && isempty(strfind(strjoin(get(H.est,'string').',' '),'SLOW')));
set(H.list,'value',5); sma_gui('select',fig);
check('SE case: p_SE.mat loaded (4 coefficients, 5th field disabled)', strcmp(get(H.pe(5),'enable'),'off') && ~isempty(get(H.pe(4),'string')));
set(H.list,'value',1); sma_gui('select',fig);

% --- validation -------------------------------------------------------------------------------------------------
set(H.tf,'string','abc'); sma_gui('run',fig);
check('bad Tf (text) rejected', ~isempty(getappdata(fig,'LastError')) && ~isempty(strfind(getappdata(fig,'LastError'),'number')));
setappdata(fig,'LastError','');
set(H.tf,'string','-1'); sma_gui('run',fig);
check('Tf <= 0 rejected', ~isempty(strfind(getappdata(fig,'LastError'),'greater than 0')));
setappdata(fig,'LastError','');
set(H.tf,'string','0.05'); set(H.pe(3),'string','1,2x'); sma_gui('run',fig);
check('bad p coefficient rejected', ~isempty(strfind(getappdata(fig,'LastError'),'p3')));
setappdata(fig,'LastError','');
sma_gui('reset',fig);
check('Reset restores default text', strcmp(get(H.pe(3),'string'), S.pstr{3}));

% --- run 1: base, Tf = 0.05, compare ON, default p, figures show+save -------------------------------------------
set(H.tf,'string','0.05'); set(H.cmp,'value',1); set(H.figs,'value',1);
sma_gui('run',fig);
res = getappdata(fig,'LastResult');
check('run base+compare finished', ~isempty(res) && isfield(res,'reference'), sprintf('(%s)', strjoin({get(H.status,'string'){:}},' ')));
check('GUI result == scripted run_case (default p, bit-identical)', isequal(res.midspan, ref.midspan));
nfig = numel(findobj(0,'type','figure'));
check('midspan window shown in addition to the start window', nfig >= 2, sprintf('(%d figures)',nfig));
check('PNG saved', exist(fullfile(root,'results','base_midspan.png'),'file')==2);
if ~isempty(shotdir), shot(fig,fullfile(shotdir,'gui_after_run.png')); end

% --- run 2: modified p ----------------------------------------------------------------------------------------------------
setappdata(fig,'LastResult',[]);
set(H.pe(4),'string',sprintf('%.8g',d0.p(4)*0.5)); sma_gui('pedit',fig);
set(H.cmp,'value',0);
if ~isempty(shotdir), shot(fig,fullfile(shotdir,'gui_screenshot.png')); end
sma_gui('run',fig);
res2 = getappdata(fig,'LastResult');
check('run with modified p finished', ~isempty(res2) && res2.cfg.p_modified);
check('modified p changes the result', ~isequal(res2.midspan, ref.midspan) && abs(res2.cfg.p(4)-d0.p(4)*0.5) <= 1e-7*abs(d0.p(4)) && isequal(res2.cfg.p([1 2 3 5]), d0.p([1 2 3 5])));
check('PNG with _pmod suffix saved (default files untouched)', exist(fullfile(root,'results','base_pmod_midspan.png'),'file')==2);
pfile_after = dir(fullfile(root,'data','p.mat')); d1 = load(fullfile(root,'data','p.mat'));
check('data/p.mat not modified', isequal(d0.p,d1.p) && pfile_before.datenum==pfile_after.datenum);

% --- run 3: linear case through the window ------------------------------------------------------------------------------
set(H.list,'value',3); sma_gui('select',fig); setappdata(fig,'LastResult',[]);
set(H.figs,'value',3);                                  % save only
sma_gui('run',fig);
res3 = getappdata(fig,'LastResult');
check('linear case runs', ~isempty(res3) && numel(res3.TT) > 10);
% --- run 4: fields plot ----------------------------------------------------------------------------------------------------------
set(H.list,'value',1); sma_gui('select',fig); set(H.tf,'string','0.03'); set(H.fld,'value',1); set(H.figs,'value',3);
setappdata(fig,'LastResult',[]); sma_gui('run',fig);
check('field plot run (save only) finished', ~isempty(getappdata(fig,'LastResult')) && exist(fullfile(root,'results','base_fields.png'),'file')==2 && exist(fullfile(root,'results','base_strain.png'),'file')==2);

% --- run 5: animation: checkboxes 'Animate beam' + 'Save GIF' ----------------------------------------------------------
check('checkboxes for the animations (beam, SMA vs linear, martensite), Save GIF and the dashboard exist', ishandle(H.anim) && ishandle(H.canim) && ishandle(H.fanim) && ishandle(H.gif) && ishandle(H.dash) && strcmp(get(H.gif,'string'),'Save GIF'));
gfile = fullfile(root,'results','base_animation.gif');  if exist(gfile,'file'), delete(gfile); end
set(H.list,'value',1); sma_gui('select',fig); set(H.tf,'string','0.03'); set(H.figs,'value',3);
set(H.gif,'value',1); set(H.anim,'value',0); setappdata(fig,'LastResult',[]);
nfig1 = numel(findobj(0,'type','figure'));
sma_gui('run',fig); r5 = getappdata(fig,'LastResult');
check('Save GIF only: results/base_animation.gif written, no animation window', ~isempty(r5) && isfield(r5,'animation') && exist(gfile,'file')==2 && numel(findobj(0,'type','figure')) == nfig1, ...
    sprintf('%d frames, %.0f kB', r5.animation.nframes, dirbytes(gfile)/1e3));
delete(gfile);  set(H.anim,'value',1); set(H.gif,'value',0); setappdata(fig,'LastResult',[]);
sma_gui('run',fig); r6 = getappdata(fig,'LastResult');
d6 = getappdata(r6.dashboard,'dash');
check('Animate beam only: loops in a dashboard tab (Beam animation), no GIF file, no extra window', ~isempty(r6) && isfield(r6,'dashboard') && strcmp(strjoin(d6.names,'|'),'Beam animation') && ~isfield(r6,'animation') && exist(gfile,'file')~=2 && numel(findobj(0,'type','figure')) == nfig1+1);
close(setdiff(findobj(0,'type','figure'), fig));
% animation speed field: used for the GIF and as the start speed of the on-screen animation
set(H.anim,'value',1); set(H.gif,'value',1); set(H.aspeed,'string','0.5'); setappdata(fig,'LastResult',[]);
sma_gui('run',fig); r6b = getappdata(fig,'LastResult');
check('Animation speed field 0.5: GIF speed and start speed of the animation tab', ~isempty(r6b) && abs(r6b.animation.speed-0.5) < 1e-12 && abs(getappdata(r6b.dashboard,'play_speed')(1)-0.5) < 1e-12);
if exist(gfile,'file'), delete(gfile); end
close(setdiff(findobj(0,'type','figure'), fig));
set(H.aspeed,'string','5'); setappdata(fig,'LastResult',[]); sma_gui('run',fig);
check('Animation speed field out of range (5): refused, nothing run', isempty(getappdata(fig,'LastResult')));
set(H.aspeed,'string','0.1'); set(H.anim,'value',0); set(H.gif,'value',0);

% --- run 6: dashboard (default) vs separate windows; comparison + martensite animations from the window ------------------------
set(H.list,'value',1); sma_gui('select',fig); set(H.tf,'string','0.07'); set(H.fld,'value',1); set(H.cmp,'value',1); set(H.figs,'value',2);   % show only
check('dashboard checkbox is on by default', get(H.dash,'value') == 1);
close(setdiff(findobj(0,'type','figure'), fig));
setappdata(fig,'LastResult',[]); sma_gui('run',fig); r7 = getappdata(fig,'LastResult');
nd = numel(findobj(0,'tag','sma_dashboard'));  nall = numel(findobj(0,'type','figure'));
check('One window (tabs): a single dashboard window instead of separate result windows', ~isempty(r7) && nd == 1 && nall == 2, sprintf('%d dashboard, %d figures in all', nd, nall));
if ~isempty(shotdir) && nd == 1
    dd = getappdata(r7.dashboard,'dash'); dashboard_select(r7.dashboard,numel(dd.names)); figure(r7.dashboard); drawnow; pause(1);
    [st,out] = system('xdotool search --name "results" | head -1'); wid = strtrim(out);
    if st==0 && ~isempty(wid)
        [~,g] = system(sprintf('xdotool getwindowgeometry --shell %s',wid));
        system(sprintf('/usr/bin/python3 %s %s %s "%s"', fullfile(fileparts(mfilename('fullpath')),'xshot.py'), getenv('DISPLAY'), fullfile(shotdir,'gui_dashboard.png'), strtrim(strrep(g,sprintf('\n'),' '))));
    end
end
close(setdiff(findobj(0,'type','figure'), fig));
set(H.dash,'value',0); setappdata(fig,'LastResult',[]); sma_gui('run',fig); r8 = getappdata(fig,'LastResult');
check('One window unticked: separate windows as before', ~isempty(r8) && numel(findobj(0,'tag','sma_dashboard')) == 0 && numel(findobj(0,'type','figure')) >= 4);
close(setdiff(findobj(0,'type','figure'), fig));
set(H.dash,'value',1); set(H.cmp,'value',0); set(H.fld,'value',0); set(H.figs,'value',3);
set(H.canim,'value',1); set(H.fanim,'value',1); set(H.gif,'value',1); set(H.anim,'value',0);
gc = fullfile(root,'results','base_animation_compare.gif'); gfl = fullfile(root,'results','base_fields_animation.gif');
if exist(gc,'file'), delete(gc); end;  if exist(gfl,'file'), delete(gfl); end
setappdata(fig,'LastResult',[]); sma_gui('run',fig); r9 = getappdata(fig,'LastResult');
d9 = getappdata(r9.dashboard,'dash');
check('SMA vs linear + martensite animations with Save GIF: both GIF files written, both in the one window', ~isempty(r9) && exist(gc,'file')==2 && exist(gfl,'file')==2 && isfield(r9,'animation_compare') && isfield(r9,'animation_fields') && isequal(d9.names,{'SMA vs linear animation','Martensite animation'}) && numel(findobj(0,'type','figure')) == 2);
if exist(gc,'file'), delete(gc); end;  if exist(gfl,'file'), delete(gfl); end
close(setdiff(findobj(0,'type','figure'), fig));
set(H.canim,'value',0); set(H.fanim,'value',0); set(H.gif,'value',0);

% --- real mouse clicks on the controls of the running animation (X11 + xdotool) -------------------------------------------------
[st,out] = system('xdotool search --name "SMA beam - start window" | head -1');
if isunix && st == 0 && ~isempty(strtrim(out)) && ~isempty(getenv('DISPLAY')) && ~isempty(shotdir)
    flag = [tempname() '.flag'];  shotf = fullfile(shotdir,'gui_animation_tab.png');
    setappdata(0,'play_log',zeros(0,3));
    setappdata(0,'sma_play_test',struct('maxtime',40,'hook',@(f,i) play_log(f,i)));
    system(sprintf('bash %s %s %s %s > /dev/null 2>&1 &', fullfile(fileparts(mfilename('fullpath')),'xdrive_dashboard.sh'), getenv('DISPLAY'), shotf, flag));
    set(H.list,'value',1); sma_gui('select',fig); set(H.tf,'string','0.05'); set(H.fld,'value',0); set(H.cmp,'value',0); set(H.figs,'value',2);
    set(H.anim,'value',1); set(H.canim,'value',1); set(H.gif,'value',0); set(H.dash,'value',1); set(H.aspeed,'string','0.1');
    setappdata(fig,'LastResult',[]);  t0 = tic;  sma_gui('run',fig);  r10 = getappdata(fig,'LastResult');  t10 = toc(t0);
    L = getappdata(0,'play_log');
    okc = ~isempty(r10) && ishandle(r10.dashboard);
    if okc
        spd = getappdata(r10.dashboard,'play_speed');  non = getappdata(r10.dashboard,'play_on');  dd10 = getappdata(r10.dashboard,'dash');
    end
    check('clicks: Stop ends the playback loop (run returns), the window stays open, all tabs paused', okc && ~any(non) && t10 < 38, sprintf('returned after %.1f s', t10));
    check('clicks: "+" twice raised the speed of the tab (0.1 -> 0.225), only that tab', okc && abs(spd(2)-0.225) < 1e-9 && abs(spd(3)-0.1) < 1e-9, sprintf('speeds %s', mat2str(spd,3)));
    check('clicks: tab buttons switch the animated tab (both animation tabs were played, hidden one paused)', all(ismember([2 3], unique(L(:,2)))) && numel(L) > 30, sprintf('%d loop passes, tabs %s', size(L,1), mat2str(unique(L(:,2))')));
    if okc
        setappdata(0,'sma_play_test',struct('maxtime',40,'hook',@(f,i) play_hook(f,i,'close',40)));
        fid = fopen(flag,'w'); fclose(fid);
        w0 = tic;  while ishandle(r10.dashboard) && toc(w0) < 20, pause(0.2); end
        check('clicks: Play restarts the loop after Stop; closing the window then ends it cleanly', ~ishandle(r10.dashboard), sprintf('%.1f s', toc(w0)));
    end
    delete(flag);
else
    fprintf('SKIP  clicks on the animation controls (needs X11 + xdotool and a screenshot folder)\n');
end
rmappdata(0,'sma_play_test');
close all
if nfail == 0, fprintf('\nALL GUI TESTS PASSED\n'); else, fprintf('\n%d GUI TEST(S) FAILED\n',nfail); end
end

function b = dirbytes(f)
d = dir(f); b = d.bytes;
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end

function shot(fig,fn)
% screenshot of the window via the X server (Linux, needs PIL); silently skipped otherwise
% bring the start window to the front: close the result windows first
for h = findobj(0,'type','figure')', if h ~= fig, close(h); end, end
figure(fig); drawnow; pause(1);
if isunix && ~isempty(getenv('DISPLAY'))
    [st,out] = system(sprintf('xdotool search --name "SMA beam - start window" | head -1'));
    wid = strtrim(out);
    if st==0 && ~isempty(wid)
        [~,g] = system(sprintf('xdotool getwindowgeometry --shell %s',wid));
        system(sprintf('/usr/bin/python3 %s %s %s "%s"', fullfile(fileparts(mfilename('fullpath')),'xshot.py'), getenv('DISPLAY'), fn, strtrim(strrep(g,sprintf('\n'),' '))));
    end
end
end
