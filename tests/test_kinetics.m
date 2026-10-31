% Tests for kinetics and rate laws

setup_test_path;
network = example_network;
s = [1 1 2 1.5 0.5 0.5]';

%% modular rate laws give finite rates
for type = {'cs','ms','ds','rp','fd'}
  network.kinetics = set_kinetics(network, type{1});
  v = network_velocities(s, network);
  assert(numel(v) == 5 && all(isfinite(v)), sprintf('kinetics type %s', type{1}));
end

%% mass-action kinetics ('standard' kinetics type)
k.type = 'standard';
for i = 1:5, k.reactions{i} = struct('type','mass-action','k',2,'km',1); end
v_expected = 2 * [1 1 2 1.5 1.5]' - [2 2 1.5 0.5 0.5]';
v = network_velocities(s, network, k);
assert(max(abs(v - v_expected)) < 1e-12);
v_split = network_velocities(s, network, k, 1);
assert(max(abs(v_split(1:2:end) - v_split(2:2:end) - v_expected)) < 1e-12);

%% kinetics types that are not implemented give an error
try
  set_kinetics(network, 'numeric');
  error('no error');
catch err
  assert(contains(err.message, 'cannot be set with set_kinetics'));
end

%% network_velocities: the argument indices selects reactions for all kinetics types
for kinetics_type = {'cs', 'ms', 'ds', 'rp', 'fd', 'mass-action', 'convenience'}
  network.kinetics = set_kinetics(network, kinetics_type{1});
  v_all    = network_velocities(ones(6,1), network);
  v_chosen = network_velocities(ones(6,1), network, network.kinetics, 0, [2 3]);
  assert(isequal(v_chosen, v_all([2 3])));
end

%% kinetics of type 'standard' give finite rates
network.kinetics = set_kinetics(network, 'standard');
assert(all(isfinite(network_velocities(s, network))));
