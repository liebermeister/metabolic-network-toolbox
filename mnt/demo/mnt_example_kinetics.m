% MNT_EXAMPLE_KINETICS Rate laws, steady states, and dynamic simulation
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

% Choose common saturable rate laws for all reactions (see 'help kinetics')

network.kinetics = set_kinetics(network, 'cs');

% Compute reaction rates for given metabolite concentrations

s = [1 1 0.5 0.5 0.1 0.1]';
v = network_velocities(s, network)

% Compute the steady state

[S, J] = network_steady_state(network, s)

% Simulate a time course of 5 time units

[t, s_t] = network_integrate(network, s, 5);
plot(t, s_t); xlabel('Time'); ylabel('Concentration'); legend(network.metabolites);
