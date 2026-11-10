function [ status ] = funforce( T , ~ , ~ )
global FF
if T > 0.1
    FF = 0*FF;
end
status = 0;
end

