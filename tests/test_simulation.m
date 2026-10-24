% Tests for simulation and steady states

setup_test_path;
network = example_network;
s0 = [1 1 0.5 0.5 0.2 0.2]';

%% network_integrate and network_integrate_par_pert agree (mass-action kinetics)
network.kinetics = set_kinetics(network, 'mass-action');
[~, s1] = network_integrate_par_pert(network, s0, 5, 0, {}, 0.1);
[~, s2] = network_integrate(network, s0, 5);
assert(max(abs(s1(:,end) - s2(:,end))) < 1e-3);

%% network_integrate: numeric kinetics give the same result as the kinetics they wrap
network.kinetics = set_kinetics(network, 'mass-action');
network_numeric = network;
network_numeric.kinetics = struct('type', 'numeric', ...
                                  'velocity_function', @(c, p, t) network_velocities(c, p.network), ...
                                  'parameters', struct('network', network));
[~, s1] = network_integrate(network, s0, 5);
[~, s2] = network_integrate(network_numeric, s0, 5);
assert(max(abs(s1(:,end) - s2(:,end))) < 1e-3);

%% network_steady_state: internal metabolites are stationary
network.kinetics = set_kinetics(network, 'cs');
[S, J] = network_steady_state(network, s0, 100);
internal = find(network.external == 0);
assert(max(abs(network.N(internal,:) * J)) < 1e-6);
