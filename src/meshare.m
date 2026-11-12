function meshare(epsyylon,xysp,xysm,xyT,t,thick,L,fig)
% MESHARE  strain and martensite-fraction fields (xi_s+, xi_s-, xi_T) over the beam at the LAST
% time step (thesis Fig. 'Simulazione_lunga'). Uses src/smooth_mean2d.m for smoothing.
% [2026 cleanup] thick / L are now arguments (default 20 / 1000 as in the original, which
% hard-coded the 30x20 mm section; the base beam is 15x10 mm -> pass 10).
if nargin < 6, thick = 20; end
if nargin < 7, L = 1000; end
if nargin < 8 || isempty(fig), fig = figure; end   % [2026] draw in the figure handed over by run_case (no second, empty window)
y = linspace(thick/2,-thick/2,10)';
for ind = length(xysp(1,:))
    epsilones = reshape(epsyylon(:,ind),length(y),[]);
    epsiloness = smooth_mean2d(epsilones,1,8); 
%     ss = reshape(s(:,ind),length(y),[]);
    xip = reshape(xysp(:,ind),length(y),[]);
    xips = smooth_mean2d(xip,1,8);  
    xim = reshape(xysm(:,ind),length(y),[]);
    xims = smooth_mean2d(xim,1,8);    
    xit = reshape(xyT(:,ind),length(y),[]);
    xits = smooth_mean2d(xit,1,8);    
    x = linspace(0,L,length(epsilones(1,:)))';
    [X,Y] = meshgrid(x,y);
    xx = linspace(0,L,400);
    yy = linspace(thick/2,-thick/2,200);
    [xx,yy]=meshgrid(xx,yy);
    epsilonesi = interp2(X,Y,epsiloness,xx,yy);
    xipi = interp2(X,Y,xips,xx,yy);
  
    ximi = interp2(X,Y,xims,xx,yy);

    xiti = interp2(X,Y,xits,xx,yy);
    figure(fig); clf
%     subplot(4,1,1)
%     mesh(xx,yy,epsilonesi)
%     xlabel('length [mm]')
%     ylabel('h [mm]')
%     colorbar
%     h = colorbar;
%     ylabel(h, '\epsilon')
%     set(h, 'Fontsize', 10)

%     drawnow
%     subplot(2,1,2)
%     mesh(x,y,ss)
%     xlabel('length [mm]')
%     ylabel('thickness [mm]')
%     colorbar
%     title('Strain field')
%     drawnow
    subplot(3,1,1)
    mesh(xx,yy,xipi)
    xlabel('length [mm]')
    ylabel('h [mm]')
    colorbar
    h = colorbar;
    ylabel(h, '\xi_s^+')
    set(h, 'Fontsize', 10)
    hTitle = title(sprintf('t = %.3g s', t(end)));   % [2026] was the fixed text 't = 0.8'
    set(hTitle, 'Fontsize', 14)
    subplot(3,1,2)
    mesh(xx,yy,ximi)
    xlabel('length [mm]')
    ylabel('h [mm]')
    colorbar
    h = colorbar;
    ylabel(h, '\xi_s^-')
    set(h, 'Fontsize', 10)
    subplot(3,1,3)
    mesh(xx,yy,xiti)
    xlabel('length [mm]')
    ylabel('h [mm]')
    colorbar
    h = colorbar;
    ylabel(h, '\xi_T')
    set(h, 'Fontsize', 10)
end
end
