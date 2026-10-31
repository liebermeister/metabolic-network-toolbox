% MNT_EXAMPLE_GRAPHICS Draw a network, with concentrations and fluxes
%
% For the network, see mnt_example_network

N = [ -1  0  0  0  0; ...
       0 -1  0  0  0; ...
       1  1 -1  0  0; ...
       0  0  1 -1 -1; ...
       0  0  0  1  0; ...
       0  0  0  0  1];
reversible   = [1 1 1 1 1]';
ind_external = [1 2 5 6]';
network      = network_construct(N, reversible, ind_external);
network.kinetics = set_kinetics(network, 'cs');

% Define the positions of metabolites and reactions

met_x = [0.1 0.1 0.3 0.7 0.9 0.9];
met_y = [0.8 0.2 0.5 0.5 0.8 0.2];
rea_x = [0.2 0.2 0.5 0.8 0.8];
rea_y = [0.65 0.35 0.5 0.65 0.35];

network = netgraph_make_graph(network, [], [], (1:11)', [met_x, rea_x; met_y, rea_y]);

% Draw the network

figure(1); netgraph_draw(network);

% Show the steady state on the network

[S, J] = network_steady_state(network, [1 1 0.5 0.5 0.1 0.1]');
figure(2); netgraph_concentrations(network, S, J, 1);
