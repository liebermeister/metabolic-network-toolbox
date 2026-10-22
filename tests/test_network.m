% Tests for network construction

setup_test_path;
network = example_network;

%% network_construct: sizes and external metabolites
assert(isequal(size(network.N), [6 5]));
assert(isequal(find(network.external(:))', [1 2 5 6]));
assert(numel(network.metabolites) == 6 && numel(network.actions) == 5);

%% network_choose: subnetwork
[subnetwork, indm, indr] = network_choose(network, [3 4], [3]);
assert(isequal(size(subnetwork.N), [numel(indm) numel(indr)]));
