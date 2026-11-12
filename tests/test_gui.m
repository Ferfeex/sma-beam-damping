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
check('slow case shows warning text', ~isempty(strfind(strjoin(get(H.est,'string').',' '),'SLOW')));
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
check('field plot run (save only) finished', ~isempty(getappdata(fig,'LastResult')) && exist(fullfile(root,'results','base_fields.png'),'file')==2);

close all
if nfail == 0, fprintf('\nALL GUI TESTS PASSED\n'); else, fprintf('\n%d GUI TEST(S) FAILED\n',nfail); end
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
