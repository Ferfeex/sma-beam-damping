function play_hook(fig, it, action, at)
% PLAY_HOOK  helper of tests/run_tests.m and tests/test_gui.m: at loop pass `at` of dashboard_play do `action`
%   ('close' the window, 'stop' = press Stop, 'pause' = press Pause of the active tab, 'tab' = switch to tab 1).
if it ~= at || ~ishandle(fig), return; end
k = getappdata(fig,'play_active');
switch action
    case 'close', close(fig);
    case 'stop',  dashboard_ctrl(fig,k,'stop');
    case 'pause', dashboard_ctrl(fig,k,'toggle');
    case 'tab',   dashboard_select(fig,1);
end
end
