% Tests for temporal_control_coefficients

setup_test_path;
network = example_network;
N = full(network.N);
rng(1);
Ec = -N' .* (1 + rand(size(N')));   % elasticities with realistic signs (stable state)
[nm, nr] = size(N);
internal = find(network.external == 0);

%% t = 0: no change of concentrations, flux changes equal the perturbation
[CS_t, CJ_t] = temporal_control_coefficients(N, Ec, network.external, 0);
assert(max(abs(CS_t(:))) < 1e-12);
assert(max(max(abs(CJ_t - eye(nr)))) < 1e-12);

%% large t: steady-state control coefficients
[CJ, CS] = control_coefficients(N, Ec, network.external);
[CS_t, CJ_t] = temporal_control_coefficients(N, Ec, network.external, 1e3);
assert(max(max(abs(CS_t - CS))) < 1e-8);
assert(max(max(abs(CJ_t - CJ))) < 1e-8);

%% agreement with a numerical integration of the linearised system
du = zeros(nr,1); du(1) = 0.1;                    % step change of reaction rate 1
N_int  = N(internal,:);
Ec_int = Ec(:,internal);
t_list = [0.5 2];
[~, x] = ode45(@(t,x) N_int * (du + Ec_int * x), [0 t_list], zeros(length(internal),1), odeset('RelTol',1e-10,'AbsTol',1e-12));
CS_t = temporal_control_coefficients(N, Ec, network.external, t_list);
for it = 1:length(t_list)
  assert(max(abs(CS_t(internal,:,it) * du - x(it+1,:)')) < 1e-6);
end

%% scaled coefficients
c = 1 + rand(nm,1); v = 1 + rand(nr,1);             % any positive values
[CS_t, CJ_t, CS_t_scaled, CJ_t_scaled] = temporal_control_coefficients(N, Ec, network.external, 1, c, v);
assert(max(max(abs(CS_t_scaled(internal,:) - diag(1./c(internal)) * CS_t(internal,:) * diag(v)))) < 1e-12);
assert(max(max(abs(CJ_t_scaled - diag(1./v) * CJ_t * diag(v)))) < 1e-12);

%% dilution: large t gives the steady-state control coefficients with dilution
mu = 0.01;
[CJ, CS] = control_coefficients(N, Ec, network.external, [], [], [], [], mu);
[CS_t, CJ_t] = temporal_control_coefficients(N, Ec, network.external, 1e4, [], [], mu);
assert(max(max(abs(CS_t - CS))) < 1e-8);
assert(max(max(abs(CJ_t - CJ))) < 1e-8);
