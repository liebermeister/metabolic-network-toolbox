function s = network_print_kinetics(network)

% NETWORK_PRINT_KINETICS Display the rate laws of a network
%
% s = network_print_kinetics(network)
%
% Display kinetic rate laws
%
% See also: set_kinetics


s = [sprintf('Network contains a kinetics of format: %s \n\n',  network.kinetics.type) ...
     kinetics_string(network.kinetics,1:length(network.actions),network.actions,network)];
