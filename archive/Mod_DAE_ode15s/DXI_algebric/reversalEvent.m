function [position, isterminal, direction] = reversalEvent(T,X)
global xisp0 xism0 xiT0
position = X(13);
isterminal = 0;
direction = 0;
if position == 0
    xisp0 = X(61:810);
    xism0 = X(811:1560);
    xiT0 = X(1561:2310);
end
end
