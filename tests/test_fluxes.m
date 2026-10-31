% Tests for flux analysis (flux balance analysis, flux sampling, thermodynamics)

setup_test_path;
network = example_network;
N = network.N;
internal = find(network.external == 0);

% flux signs, prior, and data for the tests of optimise_delta_mu_given_signs and optimise_mu_given_signs
v_signs = [1; 1; 2; 1; 1];
prior   = struct('mu0', struct('mean', zeros(6,1), 'cov', eye(6)));
data    = struct('log_c', struct('mean', zeros(6,1), 'std', ones(6,1)));
options = struct('epsilon_delta_mu', 0.1, 'max_delta_mu', 100);

c = fba_default_options(network);
c.zv    = [0 0 0 1 0]';
c.v_min = -ones(5,1);
c.v_max =  ones(5,1);

%% fba: stationary, within bounds, optimal
[v, value] = fba(network, c, 'linprog');
assert(max(abs(N(internal,:) * v)) < 1e-6);
assert(all(abs(v) <= 1 + 1e-6));
assert(abs(value - 1) < 1e-6);

%% pmf: stationary, reaches the required benefit
c.v_min(1) = 0; c.v_max(1) = 0.5;
v = pmf(network, c, 0.8);
assert(max(abs(N(internal,:) * v)) < 1e-6);
assert(abs(v(4) - 0.8) < 1e-6);

%% sample_feasible_v: samples are stationary and within bounds
[options, constraints] = stm_default_options(network);
constraints.v_min = -ones(5,1);
constraints.v_max =  ones(5,1);
V = sample_feasible_v(N, [1 2 5 6], constraints, options);
assert(max(max(abs(N(internal,:) * V))) < 1e-6);
assert(all(abs(V(:)) <= 1 + 1e-6));

%% sample_feasible_dmu: reaction Gibbs energies have the right signs
[options, constraints] = stm_default_options(network);
options.n_dmu_samples = 5;
v = [1 1 2 1 1]';
[dmu, success] = sample_feasible_dmu(N, [1 2 5 6], v, constraints, options);
assert(success == 1);
assert(all(all(dmu .* repmat(v, 1, size(dmu,2)) < 0)));

%% optimise_delta_mu_given_signs: reaction Gibbs energies have the signs required by the fluxes
delta_mu = optimise_delta_mu_given_signs(network.N, prior, data, v_signs, options);
assert(all(sign(delta_mu) == -sign(v_signs)));

%% optimise_mu_given_signs: reaction Gibbs energies have the signs required by the fluxes
constraints = struct('dmu_limit', 100);
delta_mu = optimise_mu_given_signs(network.N, prior, data, v_signs, options, constraints);
assert(all(sign(delta_mu) == -sign(v_signs)));
