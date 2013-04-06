% s = print_parameter(parameter,value)
%
%Formatted output for parameters
    
function s = print_parameter(parameter,value)

s =    sprintf('Parameter %s = %f',parameter,value);
