function [xs, ws] = strain_figure(epsyylon, XX, N, L, thick, t_end, fig)
% STRAIN_FIGURE  thesis-style picture of the beam at the LAST time step:
%   top    : vertical displacement along the beam [mm], deformed (cubic Hermite through the nodal
%            deflections w and rotations theta) against the undeformed axis, with nodes marked;
%   bottom : strain field epsilon over length x thickness (2-D colour map, same smoothing and
%            interpolation as src/meshare.m: smooth_mean2d(.,1,8) + interp2).
%   epsyylon : strain at the integration points, one column per time step (res.fields.epsyylon)
%   XX       : state history [velocities(2N-2) displacements(2N-2)] (res.XX), last row is used
%   xs, ws   : (optional outputs) abscissae and deflection of the plotted deformed shape; ws at x = L/2
%              equals res.midspan(end) (src/midspan_deflection.m, checked in tests/run_tests.m)
%   fig      : (optional) figure handle to draw in; a new figure is opened if omitted.
if nargin < 7 || isempty(fig), fig = figure; end
figure(fig); clf
Lk = L/(N-1);
u  = XX(end, 2*N-1:end).';                 % free displacement dofs at the last step
% nodal deflection w_k and rotation th_k = dw/dx (dof map of build_beam.m, verified: element dofs
% [w_a th_a w_b th_b], node k=2..N-1 -> dofs 2k-3, 2k-2; the end rotations are dofs 2N-3 (node 1)
% and 2N-2 (node N); w = 0 at both supports; same map as src/midspan_deflection.m)
w  = zeros(N,1);  th = zeros(N,1);
for k = 2:N-1, w(k) = u(2*k-3); th(k) = u(2*k-2); end
th(1) = u(2*N-3);  th(N) = u(2*N-2);
% smooth deformed shape: Hermite interpolation inside every element
xs = []; ws = [];
for k = 1:N-1
    s = linspace(0,1,25).';  if k < N-1, s = s(1:end-1); end
    h1 = 1-3*s.^2+2*s.^3;  h2 = Lk*(s-2*s.^2+s.^3);  h3 = 3*s.^2-2*s.^3;  h4 = Lk*(-s.^2+s.^3);
    xs = [xs; (k-1)*Lk + s*Lk];
    ws = [ws; h1*w(k) + h2*th(k) + h3*w(k+1) + h4*th(k+1)];
end
xn = (0:N-1)*Lk;

ax1 = subplot(2,1,1);
plot([0 L],[0 0],'k--','LineWidth',1); hold on
plot(xs,ws,'b','LineWidth',2);
plot(xn,w,'bo','MarkerSize',4,'MarkerFaceColor','b');
grid on; xlabel('length [mm]'); ylabel('vertical displacement [mm]');
title(sprintf('Deformed beam at t = %.3g s', t_end));
legend('undeformed','deformed','nodes','Location','northwest');
xlim([0 L]);
r_ = max(max(ws)-min(0,min(ws)), 1e-12);  ylim([min(0,min(ws))-0.08*r_, max(0,max(ws))+0.65*r_]);   % head-room for the legend
set(ax1,'Position',[0.12 0.64 0.80 0.28]);

% ---- strain field ----
y = linspace(thick/2,-thick/2,10)';
eps_ = reshape(epsyylon(:,end), length(y), []);
eps_s = smooth_mean2d(eps_,1,8);
x = linspace(0,L,size(eps_,2))';
[X,Y] = meshgrid(x,y);
xx = linspace(0,L,400);  yy = linspace(thick/2,-thick/2,200);
[xx,yy] = meshgrid(xx,yy);
Ei = interp2(X,Y,eps_s,xx,yy);
ax2 = subplot(2,1,2);
imagesc(xx(1,:), yy(:,1), Ei); set(gca,'YDir','normal');
colormap(jet); h = colorbar; ylabel(h,'\epsilon'); set(h,'FontSize',10);
xlabel('length [mm]'); ylabel('thickness [mm]'); title('Strain field');
set(ax2,'Position',[0.12 0.12 0.64 0.30]);
end
