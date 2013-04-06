function values = network_get_parameters(network,names)

% values = network_get_parameters(network, names)
%
% works only for 'numeric' kinetics

values = [];

switch network.kinetics.type,
  case 'numeric',
    for it = 1:length(names),
       values(it,:) = getfield(network.kinetics.parameters,names{it});
    end
  otherwise, warning('Network_get_parameters not implemented for this type of kinetics'); 
end
