function out = dashboard_ctrl(fig, k, action, arg)
% DASHBOARD_CTRL  callbacks of the control strip of animation tab k in a dashboard window (src/dashboard.m)
%   action  'toggle' (Play/Pause), 'slower', 'faster', 'slider', 'stop', 'label' (text for speed arg), 'speed' (set speed arg)
%   The speed is in simulated seconds per real second, limited to 0.02 ... 1.
out = '';
if strcmp(action,'label'), out = speed_label_(arg); return; end
if ~ishandle(fig) || ~isappdata(fig,'play_on'), return; end
h = getappdata(fig,'play_ctrl');
sp = getappdata(fig,'play_speed');
on = getappdata(fig,'play_on');
switch action
    case 'toggle'
        on(k) = ~on(k);
        setappdata(fig,'play_on',on);
        set_buttons_(h, on);
        if on(k) && ~getappdata(fig,'play_run')
            dashboard_play(fig);                  % the loop had ended (Stop / run finished): start it again from here
        end
    case 'stop'
        on(:) = false;
        setappdata(fig,'play_on',on);
        setappdata(fig,'play_run',false);
        set_buttons_(h, on);
    case 'slider'
        set_speed_(fig, k, 10^get(h(k,2),'value'), false);
    case 'faster'
        set_speed_(fig, k, sp(k)*1.5, true);
    case 'slower'
        set_speed_(fig, k, sp(k)/1.5, true);
    case 'speed'
        set_speed_(fig, k, arg, true);
    otherwise
        error('dashboard_ctrl: unknown action "%s"', action);
end
end

function set_speed_(fig, k, v, move_slider)
h = getappdata(fig,'play_ctrl');
sp = getappdata(fig,'play_speed');
v = min(1, max(0.02, v));
sp(k) = v;
setappdata(fig,'play_speed',sp);
if move_slider, set(h(k,2),'value',log10(v)); end
set(h(k,3),'string',speed_label_(v));
end

function set_buttons_(h, on)
for i = 1:size(h,1)
    if h(i,1) ~= 0 && ishandle(h(i,1))
        if on(i), set(h(i,1),'string','Pause'); else, set(h(i,1),'string','Play'); end
    end
end
end

function s = speed_label_(v)
s = sprintf('%.3gx real time (%.3gx slow motion)', v, 1/v);
end
