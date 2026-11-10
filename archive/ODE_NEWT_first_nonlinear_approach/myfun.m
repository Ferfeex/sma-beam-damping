function status = myfun(T,X,flag) 
global sigma_old differ xisp0 xism0 xiT0 sigma_i xisp xism xiT tt yy epsyylon epsilon xysp xysm xyT ind_old s
if nargin < 3 || isempty(flag)
    % VARIABLES UPDATE AND EVALUATION
    tt = [tt T];
    yy = [yy X];
    differ = [differ(:,end),sigma_i-sigma_old];
    differ = differ(:,end-1:end);
    sigma_old =sigma_i;
    ind = find(differ(:,end-1).*differ(:,end)<0);
    % WHEN DIFFER PRODUCT IS NEGATIVE, THE DIRECTION OF SIGMA HAS CHANGED
    if any(ind) ~= 0  
        % RE-CALLING THE MATERIAL CONSTANTS
        Ea = 34.9e3;
        Emp = 26e3;
        Emm = 33e3;
        EmT = 26e3;
        ssp = 0;
        sfp = 220;
        ssm = 0;
        sfm = -220;
        eps_Lp = 0.0036;
        eps_Lm = -0.0036;
        % UPDATE OF THE XI_0 CONSTANTS, ACCORDING TO MODEL IN PRESENCE OF SWITCHES
        xisp0(ind) = xysp(ind,end); % update of xi0
        xism0(ind) = xysm(ind,end);
        xiT0(ind) = xyT(ind,end);
        xisp(ind) = xysp(ind,end); % Re-assigning previous value to xi
        xism(ind) = xysm(ind,end);
        xiT(ind) = xyT(ind,end);
        % SIGMA WITH OLD EPSILON AND NEW XI
        sigma_i(ind) = (Ea + xisp(ind)*(Emp-Ea) + xism(ind)*(Emm-Ea) + xiT(ind)*(EmT-Ea)) .* (epsilon(ind) - eps_Lp*xisp(ind) - eps_Lm*xism(ind));
        differ(ind,end) = 0; % X0(IND) MUST NOT UPDATE ON NEXT ITERATION
        sigma_old = sigma_i;
    end
    if any(ind_old) ~=0
        % ALSO THOSE X0(IND_OLD) MUST NOT UPDATE
        differ(ind_old,end) = 0;
    end
    ind_old = ind; % ind_old update
    %VARIABLES UPDATE
    epsyylon = [epsyylon epsilon];
    xysp = [xysp xisp];
    xysm = [xysm xism];
    xyT = [xyT xiT];
    s = [s sigma_i];
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
