function [wm, w_node, info] = midspan_deflection(XX, N, L)
% MIDSPAN_DEFLECTION  vertical deflection of the beam at EXACTLY x = L/2, from the state history XX.
%
%   [wm, w_node, info] = midspan_deflection(XX, N, L)
%     XX      state history, one row per time step, X = [velocities(2N-2) ; displacements(2N-2)]
%     wm      deflection at x = L/2 [mm]: cubic Hermite interpolation with the shape functions of the
%             Euler-Bernoulli element (the same as in the FE model) through w and theta of the two nodes
%             that bracket L/2 (nodes 8 and 9 for N = 16). If a node sits exactly at L/2 (N odd) its
%             deflection is returned.
%     w_node  the OLD "mid-span" value: deflection of one node near the middle (node N/2+1 = 9 for N = 16,
%             x = 533 mm, which is what the original code and earlier versions of this code used)
%     info    struct: .x_node, .node (index), .col (column in XX of that node's deflection),
%             .elem, .s (position of L/2 inside the bracketing element, 0..1), .nodes (bracketing nodes)
%
% DOF MAP (checked against src/build_beam.m: Extr / Extr_tr; the element dofs are [w_a theta_a w_b theta_b]
% with theta = dw/dx - the element B matrix is built from the second derivatives of the same Hermite
% functions, curvature = w'' = Bk*[...], strain = -y*w''):
%   free displacement dofs u(1:2N-2) = XX(:, 2N-1:end)
%   node k = 2..N-1 : w_k = u(2k-3), theta_k = u(2k-2)
%   node 1          : w_1 = 0 (support), theta_1 = u(2N-3)
%   node N          : w_N = 0 (support), theta_N = u(2N-2)
Lk = L/(N-1);
u  = XX(:, 2*N-1:end);
nr = size(XX,1);
w  = zeros(nr,N);  th = zeros(nr,N);
for k = 2:N-1
    w(:,k)  = u(:,2*k-3);
    th(:,k) = u(:,2*k-2);
end
th(:,1) = u(:,2*N-3);
th(:,N) = u(:,2*N-2);

xm = L/2;
e  = min(floor(xm/Lk + 1e-9) + 1, N-1);        % element that contains L/2
s  = (xm - (e-1)*Lk)/Lk;                        % 0..1 inside that element
if s < 1e-9, s = 0; end                         % a node sits exactly at L/2
h1 = 1-3*s^2+2*s^3;   h2 = Lk*(s-2*s^2+s^3);    % Hermite shape functions at s
h3 = 3*s^2-2*s^3;     h4 = Lk*(-s^2+s^3);
wm = h1*w(:,e) + h2*th(:,e) + h3*w(:,e+1) + h4*th(:,e+1);

if mod(N,2) == 0, node = N/2+1; else, node = (N+1)/2; end   % N even: node 9 of 16 (as before); N odd: exact middle node
w_node = w(:,node);
info = struct('x_node',(node-1)*Lk, 'node',node, 'col',2*N-2+2*node-3, 'elem',e, 's',s, 'nodes',[e e+1]);
end
