% Tests for network_remove_zero_fluxes

setup_test_path;
network = example_network;
v = [1 0 1 0.5 0.5]';                          % reaction 2 carries no flux
ind_keep = v ~= 0;

%% reactions removed, metabolites kept
[v2, network2, ind_keep2] = network_remove_zero_fluxes(v, network);
assert(isequal(ind_keep2(:)', [1 0 1 1 1]));
assert(isequal(v2, v(ind_keep)));
assert(isequal(full(network2.N), full(network.N(:, ind_keep))));
assert(numel(network2.actions) == 4 && numel(network2.metabolites) == 6);

%% modular rate laws are reduced
network.kinetics = set_kinetics(network, 'cs');
[~, network3] = network_remove_zero_fluxes(v, network);
assert(numel(network3.kinetics.u) == 4 && isequal(size(network3.kinetics.KM), [4 6]));
assert(isequal(network3.kinetics.KV(:), network.kinetics.KV(ind_keep)));
assert(numel(network3.kinetics.Kcatf) == 4 && numel(network3.kinetics.Kcatr) == 4);
c = network.kinetics.c;
v_full = network_velocities(c, network);
v_red  = network_velocities(c, network3);
assert(max(abs(v_red - v_full(ind_keep))) < 1e-12);   % same rates for the remaining reactions

%% kinetics with wrong sizes or without type are removed
network.kinetics = struct('Keq', ones(3,1));
lastwarn('');
warning('off'); [~, network4] = network_remove_zero_fluxes(v, network); warning('on');
assert(~isfield(network4, 'kinetics'));
