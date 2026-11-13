% Tests for the standard concentration of thermodynamic quantities (table attribute StandardConcentration)
% (conversion functions, SBtab writers and readers, prior tables)

setup_test_path;
has_sbtab = exist('sbtab_version','file') > 0;

% R1: A + B -> C (sum of stoichiometric coefficients m = -1)
% R2: C -> D + H2O (m = 0, since water is not counted)
N = [-1 0; -1 0; 1 -1; 0 1; 0 1];
network = network_construct(N, [1;1], [1;4;5], {'A';'B';'C';'D';'H2O'}, {'R1';'R2'});
keq_mM  = [2; 3];
L       = log(1000);
keq_M   = [2000; 3];   % expected values for the standard concentration 1 M

%% conversion of equilibrium constants between 1 mM and 1 M
assert(max(abs(standard_concentration_convert(keq_mM, 'Keq', network, '1 mM', '1 M') ./ keq_M - 1)) < 1e-12);
assert(max(abs(standard_concentration_convert(keq_M, 'Keq', network, '1M', 'mM') - keq_mM)) < 1e-12);
assert(max(abs(standard_concentration_convert(log(keq_mM), 'ln Keq', network, '1 mM', '1 M') - log(keq_M))) < 1e-12);
assert(isequal(standard_concentration_convert(keq_mM, 'Keq', network, '1 mM', '1mM'), keq_mM));

%% several columns (e.g. metabolic states) are converted in the same way
K = standard_concentration_convert([keq_mM, 2*keq_mM], 'Keq', network, '1 mM', '1 M');
assert(max(max(abs(K ./ [keq_M, 2*keq_M] - 1))) < 1e-12);

%% standard reaction Gibbs free energies and standard chemical potentials agree with Keq
mu0_mM  = [1; 2; -3; 4; -5];
dmu0_mM = N' * mu0_mM;
dmu0_M  = standard_concentration_convert(dmu0_mM, 'dmu0', network, '1 mM', '1 M');
mu0_M   = standard_concentration_convert(mu0_mM, 'mu0', network, '1 mM', '1 M');
assert(max(abs(N' * mu0_M - dmu0_M)) < 1e-10);
assert(abs(mu0_M(5) - mu0_mM(5)) < 1e-12);                        % water is not converted
assert(max(abs(mu0_M(1:4) - (mu0_mM(1:4) + RT * L))) < 1e-10);    % mu0(1 M) = mu0(1 mM) + RT ln 1000
keq_from_dmu0 = exp(-dmu0_M / RT) ./ exp(-dmu0_mM / RT);
assert(max(abs(keq_from_dmu0 ./ (keq_M ./ keq_mM) - 1)) < 1e-10);

%% reaction Gibbs free energies do not change
assert(isequal(standard_concentration_convert([1; 2], 'dmu', network, '1 mM', '1 M'), [1; 2]));

%% unknown standard concentrations give an error
failed = 0; try, standard_concentration_in_mM('1 uM'); catch, failed = 1; end
assert(failed);
[c, label] = standard_concentration_in_mM('1 M');
assert(c == 1000 && strcmp(label, '1M'));

%% network state file: equilibrium constants written in 1 M are read back in 1 mM
if has_sbtab
  [~, d] = evalc('network_to_sbtab(network, struct(''keq'', keq_mM, ''dmu'', [-1; -2], ''standard_concentration'', ''1 M''))');
  assert(strcmp(d.tables.EquilibriumConstant.attributes.StandardConcentration, '1M'));
  assert(strcmp(d.tables.ReactionDeltaG.attributes.StandardConcentration, '1M'));
  values = cell_string2num(sbtab_table_get_column(d.tables.EquilibriumConstant, 'Value'));
  assert(max(abs(values ./ [2000; 3] - 1)) < 1e-6);
  filename = tempname;
  evalc('sbtab_network_state_save(network, ones(5,1), [1;1], filename, [-1; -2], struct(''verbose'',0,''standard_concentration'',''1 M''), keq_mM)');
  [~, ~, ~, ~, dmu, keq] = evalc('sbtab_network_state_load([filename ''.tsv''])');
  delete([filename '.tsv']);
  assert(max(abs(keq ./ keq_mM - 1)) < 1e-6);
  assert(isequal(dmu, [-1; -2]));
end

%% table of kinetic constants: written in 1 M, read back in 1 mM; default 1 mM is declared
if has_sbtab
  network_k = network;
  network_k.kinetics = set_kinetics(network, 'cs');
  network_k.kinetics.Keq = keq_mM;
  [~, t1] = evalc('modular_rate_law_to_sbtab(network_k, [], struct(''use_sbml_ids'', 0, ''write_concentrations'', 0, ''write_enzyme_concentrations'', 0))');
  assert(strcmp(t1.attributes.StandardConcentration, '1mM'));
  [~, t2] = evalc('modular_rate_law_to_sbtab(network_k, [], struct(''use_sbml_ids'', 0, ''write_concentrations'', 0, ''write_enzyme_concentrations'', 0, ''standard_concentration'', ''1 M''))');
  assert(strcmp(t2.attributes.StandardConcentration, '1M'));
  [~, kinetics] = evalc('sbtab_to_modular_rate_law(network, t2, struct(''use_sbml_ids'', 0))');
  assert(max(abs(kinetics.Keq ./ keq_mM - 1)) < 1e-4);
end

%% prior tables: the standard concentration is read; prior values are shifted to 1 mM
if has_sbtab
  parameter_prior = parameter_balancing_prior;
  assert(strcmp(parameter_prior.standard_concentration, '1mM'));
  assert(isequal(prior_standard_concentration_shift(parameter_prior, 'Keq', network), 0));
  parameter_prior.standard_concentration = '1M';
  shift = prior_standard_concentration_shift(parameter_prior, 'Keq', network);
  assert(max(abs(shift - [-L; 0])) < 1e-12);
  shift = prior_standard_concentration_shift(parameter_prior, 'mu0', network);
  assert(max(abs(shift - [-RT*L*ones(4,1); 0])) < 1e-10);
  assert(isequal(prior_standard_concentration_shift(parameter_prior, 'KM', network), 0));
  % missing bounds of kinetic data are filled with shifted prior bounds
  kd.Keq = struct('lower', nan(2,1), 'upper', nan(2,1), 'lower_ln', nan(2,1), 'upper_ln', nan(2,1), ...
                  'mean', nan(2,1), 'std', nan(2,1), 'median', nan(2,1), 'mean_ln', nan(2,1), 'std_ln', nan(2,1));
  kd = kinetic_data_complete(kd, parameter_prior, 1, network, 0);
  lower = eval(parameter_prior.LowerBound{parameter_prior.symbol_index.Keq});
  assert(max(abs(kd.Keq.lower ./ (lower * [1/1000; 1]) - 1)) < 1e-10);
  assert(max(abs(kd.Keq.mean_ln - (log(eval(parameter_prior.PriorMedian{parameter_prior.symbol_index.Keq})) + [-L; 0]))) < 1e-10);
end

%% kinetic data file: written in 1 M, read back in 1 mM
if has_sbtab
  parameter_prior = parameter_balancing_prior;
  kd = kinetic_data_load({'equilibrium constant'}, parameter_prior, network, [], struct('use_sbml_ids', 0, 'use_kegg_ids', 0, 'verbose', 0));
  kd.Keq.median = keq_mM; kd.Keq.mean = keq_mM; kd.Keq.std = 0.1 * keq_mM; kd.Keq.geomstd = [1.2; 1.2];
  kd.Keq.lower = keq_mM / 10; kd.Keq.upper = keq_mM * 10;
  kd.Keq.mean_ln = log(keq_mM); kd.Keq.std_ln = log([1.2; 1.2]); kd.Keq.lower_ln = log(kd.Keq.lower); kd.Keq.upper_ln = log(kd.Keq.upper);
  [t, M] = kinetic_data_to_sbtab(kd, network, [], '1 M');
  assert(strcmp(t.attributes.StandardConcentration, '1M') && strcmp(t.attributes.TableID, 'ParameterData'));
  filename = [tempname '.tsv'];
  mytable(M, 0, filename);
  kd2 = kinetic_data_load({'equilibrium constant'}, parameter_prior, network, filename, struct('use_sbml_ids', 0, 'use_kegg_ids', 0, 'verbose', 0));
  delete(filename);
  assert(max(abs(kd2.Keq.mean ./ keq_mM - 1)) < 1e-4);   % (the median is recomputed from mean and std)
  [t, ~] = kinetic_data_to_sbtab(kd, network);
  assert(strcmp(t.attributes.StandardConcentration, '1mM'));
end
