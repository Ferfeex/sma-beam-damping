function dashboard_select(fig, k)
% DASHBOARD_SELECT  show tab k of the dashboard (callback of its tab buttons)
d = getappdata(fig,'dash');
setappdata(fig,'play_active',k);              % the animation player (dashboard_play) only advances the visible tab
for i = 1:numel(d.panels)
    set(d.panels(i),'visible',ternary_(i==k,'on','off'));
    set(d.buttons(i),'value',double(i==k));
end
end
function v = ternary_(c,a,b)
if c, v = a; else, v = b; end
end
