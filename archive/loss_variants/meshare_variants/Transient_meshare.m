function meshare(epsyylon,xysp,xysm,xyT,t)
L = 1000;
y = linspace(10,-10,10)';
for ind = length(epsyylon(1,:))
    epsilones = reshape(epsyylon(:,ind),length(y),[]);
    epsiloness = smooth2a(epsilones,1,8); 
%     ss = reshape(s(:,ind),length(y),[]);
    xip = reshape(xysp(:,ind),length(y),[]);
    xips = smooth2a(xip,1,8);  
    xim = reshape(xysm(:,ind),length(y),[]);
    xims = smooth2a(xim,1,8);    
    xit = reshape(xyT(:,ind),length(y),[]);
    xits = smooth2a(xit,1,8);    
    x = linspace(0,L,length(epsilones(1,:)))';
    [X,Y] = meshgrid(x,y);
    xx = linspace(0,L,400);
    yy = linspace(10,-10,200);
    [xx,yy]=meshgrid(xx,yy);
    epsilonesi = interp2(X,Y,epsiloness,xx,yy);
    xipi = interp2(X,Y,xips,xx,yy);
  
    ximi = interp2(X,Y,xims,xx,yy);

    xiti = interp2(X,Y,xits,xx,yy);
    figure
%     subplot(4,1,1)
%     mesh(xx,yy,epsilonesi)
%     xlabel('length [mm]')
%     ylabel('thickness [mm]')
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
    xlabel('length [mm]','Fontsize', 14)
    ylabel('h [mm]','Fontsize', 14)
    colorbar
    h = colorbar;
    ylabel(h, '\xi_s^+')
    set(h, 'Fontsize', 14)
    
    hTitle = title(['t = ',num2str(t(ind)),'[s]']);
    set(hTitle, 'Fontsize', 18)
    subplot(3,1,2)
    mesh(xx,yy,ximi)
    xlabel('length [mm]','Fontsize', 14)
    ylabel('h [mm]','Fontsize', 14)
    colorbar
    h = colorbar;
    ylabel(h, '\xi_s^-')
    set(h, 'Fontsize', 14)

    subplot(3,1,3)
    mesh(xx,yy,xiti)
    xlabel('length [mm]','Fontsize', 14)
    ylabel('h [mm]','Fontsize', 14)
    colorbar
    h = colorbar;
    ylabel(h, '\xi_T')
    set(h, 'Fontsize', 14)
end
end
