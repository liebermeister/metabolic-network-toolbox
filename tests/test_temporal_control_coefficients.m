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

%% dilution with a conservation relation: warning; same coefficients as dilution of the independent metabolites only
% X (external) + B -> A + C ;  A + C -> Y (external) + B   (B + C conserved)
N2  = [-1 0; 1 -1; -1 1; 1 -1; 0 1];            % metabolites X, A, B, C, Y
ext2 = [1 0 0 0 1]';
Ec2 = [1 -0.5 1 -0.5 0; 0 1 -0.2 1 -0.3];         % substrates positive, products negative
mu = 0.05; t2 = [0.5 5 50];
lastwarn('');
CS2 = temporal_control_coefficients(N2, Ec2, ext2, t2, [], [], mu);
[~, wid] = lastwarn;
assert(strcmp(wid, 'mnt:temporal_control_coefficients:dilution_conservation'));
ii = find(~ext2); [L2, NR2] = reduce_N(N2(ii,:));
M2 = NR2 * Ec2(:,ii) * L2 - mu * eye(size(NR2,1));
for k = 1:numel(t2)
  C_red = L2 * (expm(M2 * t2(k)) - eye(size(M2))) * (M2 \ NR2);
  assert(max(max(abs(CS2(ii,:,k) - C_red))) < 1e-10);
end
lastwarn('');
temporal_control_coefficients(N2, Ec2, ext2, t2, [], [], 0);
[~, wid] = lastwarn; assert(~strcmp(wid, 'mnt:temporal_control_coefficients:dilution_conservation'));
