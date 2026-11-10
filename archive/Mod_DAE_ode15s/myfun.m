function status = myfun(T,X,flag)
global sigma_old differ xisp0 xism0 xiT0 tt yy ind_old
if nargin < 3 || isempty(flag)
    tt = [tt T];
    yy = [yy X];
    differ = [differ(:,end),X(2311:3060)-sigma_old];
    differ = differ(:,end-1:end);
    sigma_old =X(2311:3060);
    ind = find(differ(:,end-1).*differ(:,end)<=0);
    
    if any(ind)
        
        %%%%%%%%%
       Ea = 34.9e3;
       Emp = 26e3;
       Emm = 33e3;
       EmT = 26e3;
       ssp = 0;
       sfp = 188;
       ssm = 0;
       sfm = -220;
       eps_Lp = 0.016;
       eps_Lm = -0.0146;
        %%%%%%%%%%
        
        xisp0(ind) = yy(60+ind,end-1);
        xism0(ind) =  yy(810+ind,end-1);
        xiT0(ind) =  yy(1560+ind,end-1);
        X(60+ind) = xisp0;
        X(810+ind) = xism0;
        X(1560+ind) = xiT0;
        X(2310+ind) = (Ea + xisp0(ind)*(Emp-Ea) + xism0(ind)*(Emm-Ea) + xiT0(ind)*(EmT-Ea)) .* (X(3060+ind) - eps_Lp*xisp0(ind) - eps_Lm*xism0(ind));
        differ(ind,end) = 0;
    end
    if any(ind_old)
        differ(ind_old,end) = 0;
    end
    ind_old = ind;
    yy = [yy(:,1:end-1) X];
else
    switch(flag)
        case 'init'
            fprintf('start\n');
            tt = T(1);
            yy = X;
        case 'done'
            fprintf('done\n');
    end
end
status = 0;
end



