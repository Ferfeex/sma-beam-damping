function meshare(epsyylon,xysp,xysm,xyT)
L = 1000;
y = linspace(10,-10,10)';
figure
for ind = 1:50:length(xysp(1,:))
    epsilones = reshape(epsyylon(:,ind),length(y),[]);
    xip = reshape(xysp(:,ind),length(y),[]);
    xip = smooth2a(xip,1,8);  
    xim = reshape(xysm(:,ind),length(y),[]);
    xim = smooth2a(xim,1,8);    
    xit = reshape(xyT(:,ind),length(y),[]);
    xit = smooth2a(xit,1,8);    
    x = linspace(0,L,length(epsilones(1,:)))';
%     [x,y] = meshgrid(x,y);
%     xx = linspace(0,L,400);
%     yy = linspace(10,-10,200);
%     [xx,yy]=meshgrid(xx,yy);
%     epsilones = interp2(x,y,epsilones,xx,yy);
%     xip = interp2(x,y,xip,xx,yy);
  
%     xim = interp2(x,y,xim,xx,yy);

%     xit = interp2(x,y,xit,xx,yy);

    subplot(2,2,1)
    mesh(x,y,epsilones)
    xlabel('length [mm]')
    ylabel('thickness [mm]')
    colorbar
    title('Strain field')
    drawnow
    subplot(2,2,2)
    mesh(x,y,xip)
    xlabel('length [mm]')
    ylabel('thickness [mm]')
    colorbar
    title('\xi_s^+')
    drawnow
    subplot(2,2,3)
    mesh(x,y,xim)
    xlabel('length [mm]')
    ylabel('thickness [mm]')
    colorbar
    title('\xi_s^-')
    drawnow
    subplot(2,2,4)
    mesh(x,y,xit)
    xlabel('length [mm]')
    ylabel('thickness [mm]')
    colorbar
    title('\xi_T')
    drawnow
end
end
