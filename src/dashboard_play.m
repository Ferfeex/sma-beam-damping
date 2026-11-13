function info = dashboard_play(fig, varargin)
% DASHBOARD_PLAY  play the animation tabs of a dashboard window (src/dashboard.m) in a continuous loop
%   info = dashboard_play(fig, 'name', value, ...)
%
%   GNU Octave has no timer object, so the loop runs in the calling function: it BLOCKS until the window is closed or
%   Stop is pressed (run_case does this as its last step; the results are already saved at that point). In the loop
%   only the animation of the visible tab is advanced; the hidden tabs wait. Each animation wraps around at its end.
%   A frame is chosen from the elapsed real time x the speed of the tab (simulated seconds per real second,
%   0.02 ... 1, set with the slider and the -/+ buttons), so the speed does not depend on how fast the computer draws.
%   The Play/Pause button of a tab pauses only that tab. Closing the window ends the loop silently.
%
%   Options
%     'maxtime'   Inf     leave the loop after this many seconds (used by the tests)
%     'maxfps'    30      upper limit of the redraw rate
%     'hook'      []      function handle hook(fig, iteration), called once per loop pass (used by the tests)
%   info: .frames (frames drawn per tab), .wraps (completed passes per tab), .elapsed [s], .closed (window was closed)
o = struct('maxtime',Inf,'maxfps',30,'hook',[]);
for i = 1:2:numel(varargin)
    if ~isfield(o,varargin{i}), error('dashboard_play: unknown option "%s"', varargin{i}); end
    o.(varargin{i}) = varargin{i+1};
end
if isappdata(0,'sma_play_test')          % test hook: struct with .maxtime and/or .hook, used by tests/ to leave or drive the loop
    tcfg = getappdata(0,'sma_play_test');
    if isfield(tcfg,'maxtime'), o.maxtime = min(o.maxtime, tcfg.maxtime); end
    if isfield(tcfg,'hook'), o.hook = tcfg.hook; end
end
info = struct('frames',[],'wraps',[],'elapsed',0,'closed',false);
if isempty(fig) || ~ishandle(fig) || ~isappdata(fig,'play_players'), info.closed = true; return; end
if getappdata(fig,'play_run'), return; end            % a loop is already running for this window
players = getappdata(fig,'play_players');
n = numel(players);
info.frames = zeros(1,n);  info.wraps = zeros(1,n);
phase = getappdata(fig,'play_phase');
last = zeros(1,n);
setappdata(fig,'play_run',true);
t0 = tic;  tl = tic;  it = 0;
try
    while ishandle(fig) && getappdata(fig,'play_run') && toc(t0) < o.maxtime
        tw = tic;
        dt = min(toc(tl), 0.25);  tl = tic;
        k = getappdata(fig,'play_active');
        on = getappdata(fig,'play_on');
        if ~isempty(k) && k >= 1 && k <= n && ~isempty(players{k}) && on(k)
            sp = getappdata(fig,'play_speed');
            pl = players{k};
            T = pl.tf(end) - pl.tf(1);
            old = phase(k);
            phase(k) = phase(k) + dt*sp(k);
            info.wraps(k) = info.wraps(k) + floor(phase(k)/T) - floor(old/T);
            tt = pl.tf(1) + mod(phase(k), T);
            j = max(1, min(pl.nframes, round(interp1(pl.tf(:), (1:pl.nframes).', tt, 'nearest', 'extrap'))));
            if j ~= last(k)
                pl.update(j, sp(k));
                last(k) = j;
                info.frames(k) = info.frames(k) + 1;
            end
        end
        it = it + 1;
        if ~isempty(o.hook), o.hook(fig, it); end
        if ~ishandle(fig), break; end
        drawnow;
        pause(max(0.002, 1/o.maxfps - toc(tw)));
    end
catch err
    if ishandle(fig), setappdata(fig,'play_run',false); rethrow(err); end   % an error while the window is gone is ignored
end
info.elapsed = toc(t0);
info.closed = ~ishandle(fig);
if ~info.closed
    setappdata(fig,'play_run',false);
    setappdata(fig,'play_phase',phase);
end
end
