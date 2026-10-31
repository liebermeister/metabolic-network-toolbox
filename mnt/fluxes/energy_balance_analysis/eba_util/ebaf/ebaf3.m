function f = ebaf3(x,zv)

% f = ebaf3(x,zv)
%
% Linear flux benefit x(1:length(zv))'*zv (objective function used in energy balance analysis)

f = x(1:length(zv))'*zv;
