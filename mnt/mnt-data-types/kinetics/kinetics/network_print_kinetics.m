% s = network_kinetics_print(network)

function s = network_kinetics_print(network)

s = [sprintf('Network contains a kinetics of format: %s \n\n',  network.kinetics.type) ...
     kinetics_string(network.kinetics,1:length(network.actions),network.actions,network)];
