% MNT_EXAMPLE_CONTROL_ANALYSIS Elasticities, control coefficients, response coefficients
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

% Compute the steady state

[S, J] = network_steady_state(network, [1 1 0.5 0.5 0.1 0.1]');

% Elasticities and control coefficients

epsilon    = elasticities(network, S);
[C_J, C_S] = control_coefficients(network.N, epsilon, network.external);

% Check the summation theorem: this matrix must vanish

K = analyse_N(N(network.external == 0, :));
C_J * K - K

% Or compute all these quantities with a single command

R = basic_control_analysis(network, S)
