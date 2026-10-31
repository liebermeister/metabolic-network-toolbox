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

%% sample_feasible_v: lower bounds of forward fluxes are respected
N_chain = [-1 0; 1 -1; 0 1];
[~, co] = stm_default_options(N_chain);
co.v_sign = [1 1]'; co.v_min = [2 2]'; co.v_max = [10 10]';
rng(1);
v = sample_feasible_v(N_chain, [1 3]', co, struct('seed',nan,'verbose',0));
assert(all(v >= 2 - 1e-8));

%% sample_feasible_v: predefined non-zero fluxes
N_branch = [-1 0 0; 1 -1 -1; 0 1 0; 0 0 1];
[~, co] = stm_default_options(N_branch);
co.v_min = -10 * [1 1 1]'; co.v_max = 10 * [1 1 1]'; co.v_fix = [nan 3 nan]';
co.ext_signs = [-1 nan 1 1]';
rng(2);
v = sample_feasible_v(N_branch, [1 3 4]', co, struct('seed',nan,'verbose',0));
assert(abs(v(2) - 3) < 1e-8);
assert(abs(v(1) - v(2) - v(3)) < 1e-8);

%% sample_feasible_mu: rank-deficient stoichiometric matrix, all methods
N_loop = [-1 0 0 0 0; 1 -1 0 -1 0; 0 1 -1 0 0; 0 0 1 1 -1; 0 0 0 0 1];
v_loop = [2 1 1 1 2]';
[~, co] = stm_default_options(N_loop);
for method = {'sample','hit-and-run','extreme_points','centre','central'}
  [mu, feasible] = sample_feasible_mu(N_loop, [1 5]', v_loop, co, struct('seed',1), method{1}, 2);
  assert(feasible == 1);
  dmu = N_loop' * mu;
  assert(all(dmu(:) < 0) && all(abs(dmu(:)) <= co.dmu_limit + 1e-6));
end

%% sample_feasible_mu: predefined values, the other values are sampled
N_chain = [-1 0; 1 -1; 0 1];
[~, co] = stm_default_options(N_chain);
co.dmu_fix = [-3 nan]';
[mu, feasible] = sample_feasible_mu(N_chain, [1 3]', [1 1]', co, struct('seed',1), 'sample', 20);
dmu = N_chain' * mu;
assert(feasible == 1 && max(abs(dmu(1,:) + 3)) < 1e-6 && std(dmu(2,:)) > 1);

%% sample_feasible_mu: flux cycle is thermodynamically infeasible
N_cycle = [-1 0 1; 1 -1 0; 0 1 -1];
[~, co] = stm_default_options(N_cycle);
try
  [~, feasible] = sample_feasible_mu(N_cycle, [], [1 1 1]', co, struct('seed',1), 'sample', 1);
catch
  feasible = 0;
end
assert(feasible == 0);
