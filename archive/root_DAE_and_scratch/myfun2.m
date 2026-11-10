function status = myfun2(T,X,flag)
global sigma_old indp_old indm_old ex_indp ex_indm xisp_old xism_old xiT_old
if nargin < 3 || isempty(flag)
    fprintf('T = %f \n', T)
    % UPDATE OF XI_0, BEAM POINTS NOT IN TRANSFORMATION ZONE ANYMORE
    indp = find(X(2251:3000) > 0 & X(2251:3000) < 188 & X(2251:3000) >= sigma_old);
    indm = find(X(2251:3000) < 0 & X(2251:3000) > -220 & X(2251:3000) <= sigma_old);
    ex_indp = setdiff(indp_old,indp);
    ex_indm = setdiff(indm_old,indm);
    % UPDATE OF SIGMA_OLD AND IND_OLD
    if any(ex_indp)~=0 || any(ex_indm)~=0
        status = 1;
    else
        sigma_old = X(2251:3000);
        xisp_old = X(1:750);
        xism_old = X(751:1500);
        xiT_old = X(1501:2250);
        status = 0;
    end
    indp_old = indp;
    indm_old = indm;
else
    switch(flag)
        case 'init'
            fprintf('start\n');
        case 'done'
            fprintf('done\n');
    end
    status = 0;
end
end
