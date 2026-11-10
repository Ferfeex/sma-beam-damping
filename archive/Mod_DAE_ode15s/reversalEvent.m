function [position, isterminal, direction] = reversalEvent(T,X)
global xisp0 xism0 xiT0 nascar
z = double([X(15) * nascar > 0]);
position = z;
isterminal = 0;
direction = 0;
if position == 0
    xisp0 = X(61:810)
    xism0 = X(811:1560)
    xiT0 = X(1561:2310)
    T
end
end
