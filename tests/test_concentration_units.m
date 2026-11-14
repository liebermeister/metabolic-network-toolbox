% Tests for concentration units in SBtab files (conversion to mM, the unit used in MNT),
% also in combination with equilibrium constants for another standard concentration

setup_test_path;
has_sbtab = exist('sbtab_version','file') > 0;

%% unit factors
assert(concentration_unit_in_mM('M') == 1000);
assert(concentration_unit_in_mM('mM') == 1);
assert(concentration_unit_in_mM('mmol/L') == 1);
assert(abs(concentration_unit_in_mM('uM') - 0.001) < 1e-15);
assert(abs(concentration_unit_in_mM('µM') - 0.001) < 1e-15);
assert(abs(concentration_unit_in_mM('nM') - 1e-6) < 1e-18);
failed = 0; try, concentration_unit_in_mM('g/l'); catch, failed = 1; end
assert(failed == 1);

if has_sbtab

  % R1: A + B -> C (m = -1)
  N = [-1; -1; 1];
  network = network_construct(N, [1], [1;2;3], {'A';'B';'C'}, {'R1'});
  d = tempname; mkdir(d);

  %% state file in M with Keq for 1 M: concentrations -> mM, Keq -> 1 mM; driving force unchanged
  f = fullfile(d, 'state_M.tsv');
  sbtab_network_state_save(network, [0.001; 0.002; 0.0005], 1, f, [], struct('verbose',0,'standard_concentration','1 M'), 2);
  txt = fileread(f);
  % the file was written in mM (c) and 1 M (Keq); rewrite the concentrations in M
  txt = strrep(txt, 'TableName=''Metabolite concentration'' Unit=''mM''', 'TableName=''Metabolite concentration'' Unit=''M''');
  fid = fopen(f, 'w'); fprintf(fid, '%s', txt); fclose(fid);
  [~, c, ~, ~, keq] = sbtab_network_state_load(f);
  assert(max(abs(c - [1; 2; 0.5])) < 1e-12);                 % file values 0.001, 0.002, 0.0005 M, returned in mM
  assert(abs(keq - 2) < 1e-12);                               % written as Keq(1 M) = 0.002, read back as Keq(1 mM) = 2

  %% unknown concentration unit gives an error
  txt2 = strrep(txt, 'Unit=''M''', 'Unit=''g/l''');
  f2 = fullfile(d, 'state_bad.tsv'); fid = fopen(f2, 'w'); fprintf(fid, '%s', txt2); fclose(fid);
  failed = 0; try, sbtab_network_state_load(f2); catch, failed = 1; end
  assert(failed == 1);

  %% writer declares the unit mM
  f3 = fullfile(d, 'state_mM.tsv');
  sbtab_network_state_save(network, [1; 2; 3], 1, f3, [], struct('verbose',0), 2);
  assert(~isempty(strfind(fileread(f3), 'Unit=''mM''')));
  [~, c3] = sbtab_network_state_load(f3);
  assert(isequal(c3, [1; 2; 3]));

  %% kinetic data: KM values in M and uM are converted to mM
  f4 = fullfile(d, 'kinetic_data.tsv');
  fid = fopen(f4, 'w');
  fprintf(fid, '!!SBtab TableID=''KineticData'' TableType=''Quantity''\n');
  fprintf(fid, '!QuantityType\t!Reaction\t!Compound\t!Value\t!Unit\n');
  fprintf(fid, 'Michaelis constant\tR1\tA\t0.002\tM\n');
  fprintf(fid, 'Michaelis constant\tR1\tB\t30\tuM\n');
  fprintf(fid, 'Michaelis constant\tR1\tC\t0.5\tmM\n');
  fclose(fid);
  parameter_prior = parameter_balancing_prior;
  [~, kd] = evalc('kinetic_data_load({''Michaelis constant''}, parameter_prior, network, f4, struct(''use_sbml_ids'',0,''use_kegg_ids'',0,''verbose'',0))');
  km = full(kd.KM.median);
  assert(abs(km(1,1) - 2) < 1e-12);
  assert(abs(km(1,2) - 0.03) < 1e-12);
  assert(abs(km(1,3) - 0.5) < 1e-12);

  rmdir(d, 's');
end

disp('test_concentration_units: all tests passed');
