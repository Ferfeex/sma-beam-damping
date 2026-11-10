function meshare(xysp)
L = 1000;
y = linspace(10,-10,10)';
for ind = length(xysp(1,:))
    xip = reshape(xysp(:,ind),length(y),[]);
    xips = smooth2a(xip,1,8);  
    x = linspace(0,L,length(xips(1,:)))';
    [X,Y] = meshgrid(x,y);
    xx = linspace(0,L,400);
    yy = linspace(10,-10,200);
    [xx,yy]=meshgrid(xx,yy);
    xipi = interp2(X,Y,xips,xx,yy);
    figure
    mesh(xx,yy,xipi)
    xlabel('length [mm]')
    ylabel('thickness [mm]')
    colorbar
    h = colorbar;
    ylabel(h, '\xi')
    title('t = 0.8 [s]')
end
end
