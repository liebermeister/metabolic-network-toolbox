% Tests for compute_modular_elasticities and compute_modular_response_second:
% comparison with numerical derivatives of modular_velocities and with full response coefficients
% (written 2026-11-09, after the 'ds' and 'fd' branches had been repaired)

setup_test_path;

% network with a bimolecular reaction, an activator, and an inhibitor
% metabolites: 1 X (ext), 2 A, 3 B, 4 C, 5 Y (ext)
% reactions:   X -> A;   2 A <-> B + C;   B -> Y;   C -> Y
N = [-1  0  0  0;
      1 -2  0  0;
      0  1 -1  0;
      0  1  0 -1;
      0  0  1  1];
W = sparse(4,5); W(1,3) = -1; W(3,4) = 1;   % B inhibits reaction 1, C activates reaction 3
ind_ext = [1 5];
[nm, nr] = size(N);
h   = ones(nr,1);
rng(3);
c   = 0.5 + rand(nm,1);
u   = 0.5 + rand(nr,1);
KM  = sparse((0.3 + rand(nr,nm)) .* (N' ~= 0));
KA  = sparse((0.3 + rand(nr,nm)) .* (W > 0));
KI  = sparse((0.3 + rand(nr,nm)) .* (W < 0));
KV  = 1 + rand(nr,1);
Keq = exp(N' * log(c)) .* exp([2; 1.5; -1; 3]);  % driving forces of both signs (reaction 3 runs backwards)

alpha_M = sparse(nr,nm); alpha_A = sparse(nr,nm); alpha_I = sparse(nr,nm);
[ir, im] = find(KM); for k = 1:numel(ir), alpha_M(ir(k),im(k)) = KM(ir(k),im(k)) / (c(im(k)) + KM(ir(k),im(k))); end
[ir, im] = find(KA); for k = 1:numel(ir), alpha_A(ir(k),im(k)) = KA(ir(k),im(k)) / (c(im(k)) + KA(ir(k),im(k))); end
[ir, im] = find(KI); for k = 1:numel(ir), alpha_I(ir(k),im(k)) = KI(ir(k),im(k)) / (c(im(k)) + KI(ir(k),im(k))); end
A = RT * (log(Keq) - N' * log(c));

%% first- and second-order elasticities agree with numerical derivatives (all five rate laws)
for kinetic_law = {'cs','ms','ds','rp','fd'}
  kl = kinetic_law{1};
  f  = @(logc) modular_velocities(kl, N, W, ind_ext, u, exp(logc), KA, KI, KM, KV, Keq, h);
  [v, v_plus, v_minus] = f(log(c));

  E = compute_modular_elasticities(kl, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, h, v_plus, v_minus, 1);

  % first order: scaled elasticities d ln v / d ln c = (c/v) dv/dc (central differences)
  d = 1e-5; E_num = zeros(nr,nm);
  for i = 1:nm
    e = zeros(nm,1); e(i) = d;
    E_num(:,i) = (f(log(c)+e) - f(log(c)-e)) ./ (2*d) ./ v;
  end
  assert(max(abs(full(E.sc_E_c(:)) - E_num(:))) < 1e-6, sprintf('first-order elasticities wrong for %s', kl));
  un_num = diag(v) * E_num * diag(1./c);
  assert(max(abs(full(E.un_E_c(:)) - un_num(:))) < 1e-6, sprintf('unscaled elasticities wrong for %s', kl));

  % second order: d^2 ln|v| / d ln c d ln c and d^2 v / dc dc
  d = 1e-4;
  Ecc_num = zeros(nr,nm,nm); Ecc_un_num = zeros(nr,nm,nm);
  g  = @(logc) log(abs(f(logc)));
  gu = @(cc) modular_velocities(kl, N, W, ind_ext, u, cc, KA, KI, KM, KV, Keq, h);
  for i = 1:nm
    for j = 1:nm
      ei = zeros(nm,1); ei(i) = d; ej = zeros(nm,1); ej(j) = d;
      Ecc_num(:,i,j) = (g(log(c)+ei+ej) - g(log(c)+ei-ej) - g(log(c)-ei+ej) + g(log(c)-ei-ej)) / (4*d^2);
      ei = ei * c(i); ej = ej * c(j);
      Ecc_un_num(:,i,j) = (gu(c+ei+ej) - gu(c+ei-ej) - gu(c-ei+ej) + gu(c-ei-ej)) / (4*d^2*c(i)*c(j));
    end
  end
  Ecc     = double(full(E.sc_E_cc));
  Ecc_un  = double(full(E.un_E_cc));
  assert(max(abs(Ecc(:) - Ecc_num(:))) < 1e-5, sprintf('second-order elasticities wrong for %s', kl));
  assert(max(abs(Ecc_un(:) - Ecc_un_num(:))) < 1e-5, sprintf('unscaled second-order elasticities wrong for %s', kl));
end

%% compute_modular_response_second: linear and log-linear target functions
kl = 'cs';
[v, v_plus, v_minus] = modular_velocities(kl, N, W, ind_ext, u, c, KA, KI, KM, KV, Keq, h);
E   = compute_modular_elasticities(kl, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, h, v_plus, v_minus, 1);
ext = zeros(nm,1); ext(ind_ext) = 1;
[CJ, CS] = control_coefficients(N, E.un_E_c, ext);
Eu_un  = diag(v ./ u);
RSu    = CS * Eu_un;  RJu = CJ * Eu_un;
[~, ~, RSuu, RJuu] = response_coefficients_sparse(CS, E.un_E_c, Eu_un, E.un_E_cc, E.un_E_cu, E.un_E_uu);
RSuu = double(full(RSuu)); RJuu = double(full(RJuu));

% target functions, computed from the full (unscaled) response coefficients
gc = [0 1 0 2 0]'; gv = [0 1 0 0.5]';
z      = gc' * c + gv' * v;                        % linear target z = gc'*c + gv'*v
dz     = gc' * RSu + gv' * RJu;
d2z    = squeeze(tensorprod(gc', RSuu, 2, 1)) + squeeze(tensorprod(gv', RJuu, 2, 1));
R_lin  = (dz .* u') / z;
R_lin2 = diag(u) * d2z * diag(u) / z + diag(R_lin) - R_lin' * R_lin;

[R1, R2] = compute_modular_response_second(gc, gv, kl, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, h, v_plus, v_minus, CS, CJ);
assert(max(abs(full(R1(:)) - R_lin(:))) < 1e-8);
assert(max(abs(full(R2(:)) - R_lin2(:))) < 1e-8);

% log-linear target ln z = zc' * ln c + zv' * ln v
zc = [0 1 0 -1 0]'; zv = [0 2 0 0]';
RSu_sc = diag(1./c) * RSu * diag(u);  RJu_sc = diag(1./v) * RJu * diag(u);
R_log  = zc' * RSu_sc + zv' * RJu_sc;
RSuu_sc = zeros(nm,nr,nr); RJuu_sc = zeros(nr,nr,nr);
for i = 1:nm, RSuu_sc(i,:,:) = diag(u) * squeeze(RSuu(i,:,:)) * diag(u) / c(i) + diag(RSu_sc(i,:)) - RSu_sc(i,:)' * RSu_sc(i,:); end
for l = 1:nr, RJuu_sc(l,:,:) = diag(u) * squeeze(RJuu(l,:,:)) * diag(u) / v(l) + diag(RJu_sc(l,:)) - RJu_sc(l,:)' * RJu_sc(l,:); end
R_log2 = squeeze(tensorprod(zc', RSuu_sc, 2, 1)) + squeeze(tensorprod(zv', RJuu_sc, 2, 1));

[R1, R2] = compute_modular_response_second(zc, zv, kl, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, h, v_plus, v_minus, CS, CJ, 'loglinear');
assert(max(abs(full(R1(:)) - R_log(:))) < 1e-8);
assert(max(abs(full(R2(:)) - R_log2(:))) < 1e-8);
