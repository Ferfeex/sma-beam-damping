function varargout = sma_gui(action, varargin)
% SMA_GUI  interactive start window of the SMA beam simulations.
%
%   run_case            (no arguments) opens this window
%   sma_gui             the same
%   sma_gui('text')     text menu in the command window (used automatically when there is no display)
%
% In the window: choose the case (description and runtime estimate are shown), optionally edit the
% loss-factor polynomial p of THIS run (the file data/p.mat is never changed), tick "Compare with
% linear beam", optionally type a final time Tf, choose whether the figures are shown and/or saved
% (results/*.png), whether the plots go into ONE window with tabs (default) or separate windows, which
% animations to make (beam, SMA vs linear beam, martensite fractions; "Save GIF" writes the GIF files
% in results/), and press Run. Built only from classic uicontrols, so it works in GNU Octave (Qt
% toolkit, e.g. the Octave GUI on Windows) and in MATLAB.
%
% Returns the figure handle (window mode). The other actions ('run', 'select', ...) are the window's
% own callbacks; they are also used by the automatic test.
setup_paths();
if nargin == 0, action = 'open'; end
switch lower(action)
    case 'open'
        if has_display()
            fig = open_window();
            if nargout > 0, varargout{1} = fig; end
        else
            fprintf('(no graphical display available - using the text menu)\n');
            text_menu();
            if nargout > 0, varargout{1} = []; end
        end
    case 'text'
        text_menu();
        if nargout > 0, varargout{1} = []; end
    case 'select', cb_select(varargin{1});
    case 'pedit',  cb_pedit(varargin{1});
    case 'tfedit', cb_tfedit(varargin{1});
    case 'reset',  cb_reset(varargin{1});
    case 'run',    cb_run(varargin{1});
    case 'close',  close(varargin{1});
    otherwise
        error('sma_gui: unknown action "%s"', action);
end
end

% ======================================================================================
function ok = has_display()
ok = true;
if exist('OCTAVE_VERSION','builtin') ~= 0
    if isunix && ~ismac && isempty(getenv('DISPLAY')) && isempty(getenv('WAYLAND_DISPLAY'))
        ok = false; return
    end
    try
        if ~any(strcmp(available_graphics_toolkits(),'qt')), ok = false; return, end
        if ~strcmp(graphics_toolkit(),'qt'), graphics_toolkit('qt'); end
    catch
        ok = false;
    end
else
    try, ok = usejava('awt'); catch, ok = false; end
end
end

% ======================================================================================
function fig = open_window()
C = case_catalog();
S.C = C;
S.nP = 5;
S.figmodes = {'both','show','save','none'};
fig = figure('name','SMA beam - start window','numbertitle','off','menubar','none','toolbar','none', ...
    'units','pixels','position',[120 80 980 640],'resize','on','tag','sma_gui');
bg = get(fig,'color');
H.fig = fig;
mk = @(style,pos,str,varargin) uicontrol(fig,'style',style,'units','normalized','position',pos, ...
    'string',str,'backgroundcolor',bg,varargin{:});

mk('text',[0.02 0.935 0.96 0.05],'SMA beam damping - choose a simulation','fontsize',14,'fontweight','bold', ...
    'horizontalalignment','left');

% ---- left column: case, description, runtime, options --------------------------------------
mk('text',[0.02 0.895 0.45 0.035],'1. Case','fontweight','bold','horizontalalignment','left');
H.list = uicontrol(fig,'style','listbox','units','normalized','position',[0.02 0.64 0.45 0.255], ...
    'string',{C.label},'value',1,'min',1,'max',1,'backgroundcolor',[1 1 1], ...
    'callback',@(h,e) sma_gui('select',fig));
H.desc = uicontrol(fig,'style','text','units','normalized','position',[0.02 0.485 0.45 0.15], ...
    'string','','backgroundcolor',[0.97 0.97 0.92],'horizontalalignment','left');
H.est  = uicontrol(fig,'style','text','units','normalized','position',[0.02 0.405 0.45 0.07], ...
    'string','','backgroundcolor',bg,'horizontalalignment','left','fontweight','bold');
mk('text',[0.02 0.365 0.45 0.035],'2. Options','fontweight','bold','horizontalalignment','left');
mk('text',[0.02 0.325 0.20 0.04],'Final time Tf [s]:','horizontalalignment','left');
H.tf = uicontrol(fig,'style','edit','units','normalized','position',[0.22 0.33 0.10 0.045], ...
    'string','','backgroundcolor',[1 1 1],'callback',@(h,e) sma_gui('tfedit',fig));
H.tfhint = mk('text',[0.33 0.325 0.14 0.04],'(blank = default)','horizontalalignment','left');
H.cmp = mk('checkbox',[0.02 0.285 0.45 0.04],'Compare with linear beam (overlay)','value',0, ...
    'horizontalalignment','left');
H.fld = mk('checkbox',[0.02 0.248 0.45 0.04],'Plot fields (strain / martensite, incl. martensite vs time)','value',0, ...
    'horizontalalignment','left');
mk('text',[0.02 0.208 0.10 0.04],'Animate:','horizontalalignment','left');
H.anim = mk('checkbox',[0.11 0.212 0.12 0.04],'beam','value',0, ...
    'horizontalalignment','left');
H.canim = mk('checkbox',[0.22 0.212 0.14 0.04],'SMA vs linear','value',0, ...
    'horizontalalignment','left');
H.fanim = mk('checkbox',[0.36 0.212 0.12 0.04],'martensite','value',0, ...
    'horizontalalignment','left');
H.gif = mk('checkbox',[0.02 0.172 0.22 0.04],'Save GIF','value',0, ...
    'horizontalalignment','left');
H.dash = mk('checkbox',[0.25 0.172 0.22 0.04],'One window (tabs)','value',1, ...
    'horizontalalignment','left');
mk('text',[0.02 0.12 0.14 0.04],'Figures:','horizontalalignment','left');
H.figs = uicontrol(fig,'style','popupmenu','units','normalized','position',[0.14 0.125 0.33 0.045], ...
    'string',{'Show windows and save PNG','Show windows only','Save PNG only (no windows)','No figures'}, ...
    'value',1,'backgroundcolor',[1 1 1]);

% ---- right column: p polynomial, preview, run ---------------------------------------------------
H.ptitle = mk('text',[0.52 0.895 0.46 0.035],'3. Loss-factor polynomial p','fontweight','bold','horizontalalignment','left');
H.pform = mk('text',[0.52 0.855 0.46 0.04],'eta(x) = p1 x^4 + p2 x^3 + p3 x^2 + p4 x + p5','horizontalalignment','left');
H.pe = zeros(1,S.nP); H.pl = zeros(1,S.nP);
for k = 1:S.nP
    x0 = 0.52 + (k-1)*0.092;
    H.pl(k) = mk('text',[x0 0.815 0.085 0.03],sprintf('p%d',k),'horizontalalignment','center');
    H.pe(k) = uicontrol(fig,'style','edit','units','normalized','position',[x0 0.77 0.087 0.045], ...
        'string','','backgroundcolor',[1 1 1],'callback',@(h,e) sma_gui('pedit',fig));
end
H.reset = uicontrol(fig,'style','pushbutton','units','normalized','position',[0.52 0.71 0.20 0.05], ...
    'string','Reset to default','callback',@(h,e) sma_gui('reset',fig));
H.pnote = mk('text',[0.73 0.69 0.25 0.08],'','horizontalalignment','left');
H.ax = axes('parent',fig,'units','normalized','position',[0.58 0.30 0.39 0.34],'box','on');
H.run = uicontrol(fig,'style','pushbutton','units','normalized','position',[0.52 0.115 0.27 0.09], ...
    'string','Run','fontsize',13,'fontweight','bold','callback',@(h,e) sma_gui('run',fig));
H.close = uicontrol(fig,'style','pushbutton','units','normalized','position',[0.81 0.115 0.16 0.09], ...
    'string','Close','callback',@(h,e) sma_gui('close',fig));
H.status = uicontrol(fig,'style','text','units','normalized','position',[0.02 0.005 0.96 0.065], ...
    'string','Ready. Choose a case and press Run.','backgroundcolor',bg,'horizontalalignment','left');

S.H = H; S.pdef = []; S.pstr = {};
setappdata(fig,'S',S);
cb_select(fig);
end

% ======================================================================================
function cb_select(fig)
S = getappdata(fig,'S'); H = S.H;
k = get(H.list,'value'); C = S.C(k);
set(H.desc,'string',wrap({strjoin(C.text(:).',' ')},58));
set(H.tf,'string','');
set(H.tfhint,'string',sprintf('(blank = %g s)',C.Tf_def));
lin = strcmp(C.model,'linear');
set(H.cmp,'value',0,'enable',onoff(~lin));
set(H.fld,'value',0,'enable',onoff(strcmp(C.model,'two_phase')));
set(H.anim,'value',0); set(H.gif,'value',0);
set(H.canim,'value',0,'enable',onoff(~lin)); set(H.fanim,'value',0,'enable',onoff(~lin));
if lin
    S.pdef = []; S.pstr = repmat({''},1,S.nP);
    for i = 1:S.nP, set(H.pe(i),'string','','enable','off'); end
    set(H.reset,'enable','off');
    set(H.pnote,'string',wrap({'Not used: this is a linear beam (no loss factor).'},30));
    set(H.pform,'string','');
else
    root = fileparts(mfilename('fullpath'));
    d = load(fullfile(root,'data',C.p_file));
    p = d.p(:).';
    S.pdef = p; n = numel(p);
    S.pstr = cell(1,S.nP);
    for i = 1:S.nP
        if i <= n
            S.pstr{i} = sprintf('%.6g',p(i));
            set(H.pe(i),'string',S.pstr{i},'enable','on');
            set(H.pl(i),'visible','on');
        else
            S.pstr{i} = '';
            set(H.pe(i),'string','','enable','off');
        end
    end
    set(H.reset,'enable','on');
    if strcmp(C.model,'one_phase_SE')
        set(H.pform,'string','eta(x) = p1 x^3 + p2 x^2 + p3 x + p4,  x = stress [MPa]');
    else
        set(H.pform,'string','eta(x) = p1 x^4 + p2 x^3 + p3 x^2 + p4 x + p5,  x = strain (0.01 = 1 %)');
    end
    set(H.pnote,'string',wrap({sprintf('Default: data/%s. Edited values are used for this run only.',C.p_file)},30));
end
setappdata(fig,'S',S);
update_estimate(fig);
update_preview(fig);
end

function cb_pedit(fig)
update_preview(fig);
end

function cb_tfedit(fig)
update_estimate(fig);
end

function cb_reset(fig)
S = getappdata(fig,'S');
for i = 1:S.nP
    if ~isempty(S.pstr{i}), set(S.H.pe(i),'string',S.pstr{i}); end
end
update_preview(fig);
end

% ======================================================================================
function update_estimate(fig)
S = getappdata(fig,'S'); H = S.H;
C = S.C(get(H.list,'value'));
[Tf,err] = read_tf(H);
if ~isempty(err), Tf = []; end
if isempty(Tf), Tf = C.Tf_def; end
est = C.t_ref * Tf / C.Tf_def;
if get(H.cmp,'value'), est = est + 2; end
str = sprintf('Estimated runtime (Octave, 1 core): about %s', fmt_time(est));
if est > 20*60
    str = sprintf('WARNING - VERY SLOW: about %s in Octave. Use a short Tf (e.g. 0.02) to test.', fmt_time(est));
    set(H.est,'foregroundcolor',[0.8 0 0]);
else
    set(H.est,'foregroundcolor',[0 0 0]);
end
set(H.est,'string',wrap({str},56));
end

function update_preview(fig)
S = getappdata(fig,'S'); H = S.H;
C = S.C(get(H.list,'value'));
ax = H.ax;
cla(ax);
if isempty(S.pdef)
    set(ax,'visible','off');
    return
end
set(ax,'visible','on');
if strcmp(C.model,'one_phase_SE')
    x = linspace(0,150,200); sc = 1; xl = 'stress x [MPa]';
else
    x = linspace(0,0.025,200); sc = 100; xl = 'strain x [%]';
end
[p,err] = read_p(S);
set(ax,'nextplot','add');
plot(ax,x*sc,polyval(S.pdef,x),'b--','linewidth',1.5);
if isempty(err)
    plot(ax,x*sc,polyval(p,x),'r-','linewidth',1.5);
    plot(ax,C.x_ref*sc,polyval(p,C.x_ref),'ro','markersize',6);
    t = 'loss factor \eta(x): dashed = default, red = current';
else
    t = ['invalid input: ' err];
end
set(ax,'nextplot','replace');
xlabel(ax,xl); ylabel(ax,'\eta'); title(ax,t,'fontsize',9); grid(ax,'on');
end

% ======================================================================================
function [p,err] = read_p(S)
p = []; err = '';
if isempty(S.pdef), return, end
n = numel(S.pdef); p = zeros(1,n);
for i = 1:n
    str = strtrim(get(S.H.pe(i),'string'));
    if strcmp(str,S.pstr{i})              % untouched field: use the exact stored value
        p(i) = S.pdef(i); continue
    end
    v = str2double(strrep(str,',','.'));
    if ~isscalar(v) || ~isfinite(v)
        err = sprintf('p%d is not a number', i); p = []; return
    end
    p(i) = v;
end
end

function [Tf,err] = read_tf(H)
Tf = []; err = '';
s = strtrim(get(H.tf,'string'));
if isempty(s), return, end
v = str2double(strrep(s,',','.'));
if ~isscalar(v) || ~isfinite(v)
    err = 'The final time Tf must be a number (or blank for the default).';
elseif v <= 0
    err = 'The final time Tf must be greater than 0.';
else
    Tf = v;
end
end

function cb_run(fig)
S = getappdata(fig,'S'); H = S.H;
C = S.C(get(H.list,'value'));
silent = isappdata(fig,'Silent') && getappdata(fig,'Silent');
% ---- validate -------------------------------------------------------------------------------
[Tf,err] = read_tf(H);
if ~isempty(err), show_error(fig,err,silent); return, end
[p,err] = read_p(S);
if ~isempty(err)
    show_error(fig,['Loss-factor polynomial: ' err '. Fix it or press "Reset to default".'],silent); return
end
pmod = ~isempty(p) && ~isequal(p,S.pdef);
est = C.t_ref; if ~isempty(Tf), est = est*Tf/C.Tf_def; end
if est > 20*60 && ~silent
    a = questdlg(sprintf(['This run will take about %s in Octave (hours). Your Octave stays busy until it ' ...
        'finishes (Ctrl+C in the command window stops it).\n\nRun it anyway?'], fmt_time(est)), ...
        'Slow simulation','Run anyway','Cancel','Cancel');
    if ~strcmp(a,'Run anyway'), set(H.status,'string','Cancelled.'); return, end
end
args = {'compare_hle', logical(get(H.cmp,'value')) && strcmp(get(H.cmp,'enable'),'on'), ...
        'plot_fields', logical(get(H.fld,'value')) && strcmp(get(H.fld,'enable'),'on'), ...
        'figures', S.figmodes{get(H.figs,'value')}};
anim = logical(get(H.anim,'value'));  sgif = logical(get(H.gif,'value'));
canim = logical(get(H.canim,'value')) && strcmp(get(H.canim,'enable'),'on');
fanim = logical(get(H.fanim,'value')) && strcmp(get(H.fanim,'enable'),'on');
if anim && sgif, args = [args {'animate','both'}]; elseif sgif && ~canim && ~fanim, args = [args {'animate','gif'}]; elseif anim, args = [args {'animate',true}]; end
if canim, args = [args {'animate_compare', ternary(sgif,'both',true)}]; end      % 'Save GIF' = window AND GIF
if fanim, args = [args {'animate_fields', ternary(sgif,'both',true)}]; end
args = [args {'dashboard', logical(get(H.dash,'value'))}];
if ~isempty(Tf), args = [args {'Tf',Tf}]; end
if pmod, args = [args {'p',p}]; end
% ---- run ----------------------------------------------------------------------------------------
set(H.status,'string',sprintf('RUNNING "%s" ... the window is busy until it finishes; progress is printed in the command window.',C.name));
set(H.run,'enable','off'); drawnow;
ok = true;
try
    res = run_case(C.name,args{:});
catch e
    ok = false; msg = e.message;
end
if ~ishandle(fig), return, end            % window was closed meanwhile
set(H.run,'enable','on');
if ok
    k = min(6,size(res.peaks,1));
    pk = sprintf('%.2f  ',res.peaks(1:k,2));
    note = ''; if pmod, note = ' (modified p)'; end
    set(H.status,'string',wrap({sprintf('Done%s in %.1f s. Mid-span peaks [mm]: %s', note, res.runtime_s, pk)},110));
    setappdata(fig,'LastResult',res);
    drawnow;
else
    set(H.status,'string',wrap({['Error: ' msg]},110));
    show_error(fig,['The simulation failed: ' msg],silent);
end
end

function show_error(fig,msg,silent)
setappdata(fig,'LastError',msg);
S = getappdata(fig,'S');
set(S.H.status,'string',wrap({['Error: ' msg]},110));
if ~silent, errordlg(msg,'SMA beam'); end
end

% ======================================================================================
function text_menu()
C = case_catalog();
fprintf('\n=== SMA beam simulations (text menu) ===\n');
for k = 1:numel(C)
    fprintf('  %d) %-26s ~%s%s\n', k, C(k).name, fmt_time(C(k).t_ref), ternary(C(k).slow,'   ** VERY SLOW **',''));
    fprintf('       %s\n', C(k).label);
end
k = str2double(input(sprintf('Case number (1-%d, empty = quit): ',numel(C)),'s'));
if isnan(k) || k < 1 || k > numel(C) || k ~= round(k), fprintf('Cancelled.\n'); return, end
c = C(k);
fprintf('%s\n', strjoin(c.text.', sprintf('\n')));
args = {};
s = strtrim(input(sprintf('Final time Tf in s (empty = %g): ',c.Tf_def),'s'));
if ~isempty(s)
    Tf = str2double(s);
    if ~isfinite(Tf) || Tf <= 0, fprintf('Tf must be a number > 0. Cancelled.\n'); return, end
    args = [args {'Tf',Tf}];
else
    Tf = c.Tf_def;
end
est = c.t_ref*Tf/c.Tf_def;
if est > 20*60
    fprintf('WARNING: about %s expected (hours). Use a short Tf to test.\n', fmt_time(est));
    if ~strcmpi(strtrim(input('Continue? (y/n) ','s')),'y'), fprintf('Cancelled.\n'); return, end
end
if ~strcmp(c.model,'linear')
    root = fileparts(mfilename('fullpath'));
    d = load(fullfile(root,'data',c.p_file)); p = d.p(:).';
    fprintf('Loss-factor polynomial p (default): %s\n', sprintf('%.6g ',p));
    if strcmpi(strtrim(input('Change p for this run? (y/n) ','s')),'y')
        s = strtrim(input(sprintf('Enter the %d coefficients separated by spaces: ',numel(p)),'s'));
        q = sscanf(strrep(s,',','.'),'%f').';
        if numel(q) ~= numel(p) || any(~isfinite(q)), fprintf('Need exactly %d numbers. Cancelled.\n',numel(p)); return, end
        args = [args {'p',q}];
    end
    if strcmpi(strtrim(input('Compare with the linear beam? (y/n) ','s')),'y'), args = [args {'compare_hle',true}]; end
    if strcmp(c.model,'two_phase') && strcmpi(strtrim(input('Plot fields (strain / martensite)? (y/n) ','s')),'y')
        args = [args {'plot_fields',true}];
    end
end
if has_display(), fm = 'both'; else, fm = 'save'; end
run_case(c.name,args{:},'figures',fm);
end

% ======================================================================================
function s = onoff(tf)
if tf, s = 'on'; else, s = 'off'; end
end
function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end
function s = fmt_time(t)
if t < 90, s = sprintf('%.0f s',t);
elseif t < 5400, s = sprintf('%.0f min',t/60);
else, s = sprintf('%.1f h',t/3600);
end
end
function out = wrap(lines,w)
% word-wrap a cell array of lines to about w characters (cell array of strings)
out = {};
for i = 1:numel(lines)
    words = strsplit(strtrim(lines{i}),' ');
    cur = '';
    for j = 1:numel(words)
        if isempty(cur), cur = words{j};
        elseif numel(cur)+1+numel(words{j}) <= w, cur = [cur ' ' words{j}];
        else, out{end+1,1} = cur; cur = words{j};
        end
    end
    out{end+1,1} = cur;
end
end
