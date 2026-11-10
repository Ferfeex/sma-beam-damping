function status = myfunORIGINAL(T,X,flag)
global tt yy 
if nargin < 3 || isempty(flag)
    tt = [tt T];
    yy = [yy X];
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
