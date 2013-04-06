% s = network_print_kinetics(network)

function s = network_print_kinetics(network)

s = [sprintf('Network contains a kinetics of format: %s \n\n',  network.kinetics.type) ...
     kinetics_string(network.kinetics,1:length(network.actions),network.actions,network)];
