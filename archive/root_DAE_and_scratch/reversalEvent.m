function [position, isterminal, direction] = reversalEvent(T,X)
global xi_sp0 xi_sm0 xi_T0
position = X(13);
isterminal = 0;
direction = 0;
    xi_sp0 = X(61:810);
    xi_sm0 = X(811:1560);
    xi_T0 = X(1561:2310);
    T
end
