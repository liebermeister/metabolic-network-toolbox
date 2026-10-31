function s = cs_print(kinetics,n,actions,network)

% s = cs_print(kinetics,n,actions,network)
%
% Print enzyme levels and metabolite levels of a kinetics structure into a string

s = [sprintf('u: %f\n',kinetics.u); ...
     sprintf('c: %f\n',kinetics.c)];
