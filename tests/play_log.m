function play_log(fig, it)
% PLAY_LOG  helper of tests/test_gui.m: hook of dashboard_play that records the active tab of every loop pass
L = getappdata(0,'play_log');
if ishandle(fig), L(end+1,:) = [it getappdata(fig,'play_active') now*86400]; end   %#ok<AGROW>
setappdata(0,'play_log',L);
end
