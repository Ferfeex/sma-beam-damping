function status = myfun(T,X,flag)
global sigma_old
if nargin < 3 || isempty(flag)
    sigma_old = X(2311:3060);
else
    switch(flag)
        case 'init'
            fprintf('start\n');
            sigma_old = X(2311:3060);
        case 'done'
            fprintf('done\n');
    end
end
status = 0;
end
