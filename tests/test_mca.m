% Tests for metabolic control analysis

setup_test_path;
network = example_network;
N = full(network.N);
rng(1);
Ec = -N' .* (1 + rand(size(N')));   % elasticities with realistic signs

internal = find(network.external == 0);
K  = null(N(internal,:));           % stationary flux modes

%% control_coefficients: summation and connectivity theorems (unscaled)
[CJ, CS] = control_coefficients(N, Ec, network.external);
assert(max(max(abs(CJ * K - K))) < 1e-10);
assert(max(max(abs(CS * K))) < 1e-10);
assert(max(max(abs(CJ * Ec(:,internal)))) < 1e-10);

%% control_coefficients: tiny dilution rate gives (almost) the same result
[CJ,  CS ] = control_coefficients(N, Ec, network.external);
[CJ0, CS0] = control_coefficients(N, Ec, network.external, [], [], [], [], 1e-9);
assert(max(abs(CJ(:) - CJ0(:))) < 1e-6 && max(abs(CS(:) - CS0(:))) < 1e-6);

%% norm_response_coefficients: reactions not in used_J get zero rows
RS = rand(6,3); RJ = rand(5,3);
[~, RJ_norm] = norm_response_coefficients(ones(6,1), 2*ones(5,1), ones(3,1), RS, RJ, [], [], [1 1 0 1 1]');
assert(all(RJ_norm(3,:) == 0));
assert(max(abs(RJ_norm(1,:) - RJ(1,:)/2)) < 1e-12);

%% internal_elasticities: metabolite that is substrate and inhibitor of the same reaction
regulated_network = network;
regulated_network.regulation_matrix = sparse(5,6);
regulated_network.regulation_matrix(3,3) = -1;
s_test = 2 * ones(6,1);
for kinetics_type = {'cs', 'ms', 'mass-action', 'convenience'}
  regulated_network.kinetics = set_kinetics(regulated_network, kinetics_type{1});
  epsilon = full(internal_elasticities(regulated_network, s_test));
  s_pert = s_test; s_pert(3) = s_test(3) * 1.000001;
  v      = network_velocities(s_test, regulated_network);
  v_pert = network_velocities(s_pert, regulated_network);
  assert(norm(epsilon(:,1) - (v_pert-v)/(s_pert(3)-s_test(3))) < 1e-3);
end

%% control_coefficients: unstable state gives an error
Ec_unstable = -Ec;
failed = 0;
try
  control_coefficients(N, Ec_unstable, network.external);
catch
  failed = 1;
end
assert(failed == 1);

%% control_coefficients: on_unstable = 'nan' and output 'stable'
[CJ, CS, ~, ~, ~, ~, ~, stable] = control_coefficients(N, -Ec, network.external, [], [], [], [], [], 'nan');
assert(stable == 0 && all(isnan(CJ(:))) && all(isnan(CS(:))));
assert(isequal(size(CS), size(N)) && isequal(size(CJ), size(N,2) * [1 1]));
[~, ~, ~, ~, ~, ~, ~, stable] = control_coefficients(N, Ec, network.external);
assert(stable == 1);

%% control_coefficients: eigenvalues that are positive by numerical error only give no error
N1  = [-1 0; 1 -1; 0 1];
Ec1 = [0 1e-12 0; 0 0 0];
[~, ~, ~, ~, M, M_adj, ~, stable] = control_coefficients(N1, Ec1, [1 0 1]');
assert(stable == 1 && max(real(eig(M))) > 0 && max(real(eig(M_adj))) < 0);
