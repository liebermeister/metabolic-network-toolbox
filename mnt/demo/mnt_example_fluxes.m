% MNT_EXAMPLE_FLUXES Flux balance analysis and thermodynamic feasibility
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

% Flux balance analysis: maximise the flux through reaction 4

fba_constraints       = fba_default_options(network);
fba_constraints.zv    = [0 0 0 1 0]';
fba_constraints.v_min = -ones(5,1);
fba_constraints.v_max =  ones(5,1);

[v, value] = fba(network, fba_constraints)

% Among all flux distributions with the same benefit, find the one
% with the smallest sum of absolute fluxes

v_minimal = pmf(network, fba_constraints, value)

% Test whether the fluxes are thermodynamically feasible

feasible = eba_feasible_lp(v_minimal, network.N)
