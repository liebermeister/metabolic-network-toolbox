function network = network_check_regulation_matrix(network)

% network = network_check_regulation_matrix(network)

epsilon = elasticities(network,ones(size(network.metabolites)));

network.regulation_matrix = double( (epsilon~=0) .* (network.N'==0) );

