% MNT_EXAMPLE_NETWORK Build a network, inspect it, and change it
%
% The network used in all examples has six metabolites and five reactions.
% S1, S2, S5, and S6 are external metabolites:
%
%   S1             S5
%      \         /
%        S3 - S4
%      /         \
%   S2             S6

% Build the network from its stoichiometric matrix, the reversibility
% of the reactions, and the indices of external metabolites

N = [ -1  0  0  0  0; ...
       0 -1  0  0  0; ...
       1  1 -1  0  0; ...
       0  0  1 -1 -1; ...
       0  0  0  1  0; ...
       0  0  0  0  1];
reversible   = [1 1 1 1 1]';
ind_external = [1 2 5 6]';
network      = network_construct(N, reversible, ind_external);

% Show the reaction formulae

network_print_formulae(network)

% Or build the network directly from reaction formulae

network_2 = network_build_from_sum_formulae_list({'A <=> B'; 'B + C <=> D'});
network_print_formulae(network_2)

% Look up metabolites, and the reactions they take part in

ind_metabolite = network_find_metabolites(network, 'S3')
ind_reactions  = network_find_reactions(network, 'S3')

% Structural properties: stationary flux modes (columns of K) and link matrix L

[K, L] = network_analyse(network)

% Choose a part of the network: reactions 3 and 4, with their metabolites

subnetwork = network_choose(network, [], [3 4]);
network_print_formulae(subnetwork)
