% Tests for the conversion between SBtab and network structures (Metabolic Network Toolbox)

setup_test_path;
has_mnt = exist('network_construct','file') > 0;
if has_mnt
  N = [-1 0; 1 -1; 0 1];
  network = network_construct(N, [1;1], [1;3], {'A';'B';'C'}, {'R1';'R2'});
  network.kinetics = set_kinetics(network, 'cs');
end

%% network -> SBtab -> network
if has_mnt
  [~, d] = evalc('network_to_sbtab(network, struct(''c'',[1;0;2],''v'',[1;1]))');
  assert(all(ismember({'Reaction','Compound','MetaboliteConcentration','Flux'}, fieldnames(d.tables))));
  assert(~isfield(d.tables, 'Parameter'));
  assert(isequal(sbtab_table_get_column(d.tables.Compound,'IsConstant'), {'True';'False';'True'}));
  [~, network2] = evalc('sbtab_to_network(d)');
  assert(isequal(network2.metabolites, network.metabolites));
  assert(isequal(network2.actions, network.actions));
  assert(isequal(full(network2.N), N));
  assert(isequal(network2.external(:), [1;0;1]));
end

%% kinetic parameters are written only if requested
if has_mnt
  [~, d] = evalc('network_to_sbtab(network, struct(''modular_rate_law_table'',1))');
  assert(isfield(d.tables, 'Parameter'));
end

%% save and load network state
if has_mnt
  filename = tempname;
  evalc('sbtab_network_state_save(network, [1;0;2], [1;1], filename, [], struct(''verbose'',0))');
  [~, network3, c, v, dmu] = evalc('sbtab_network_state_load([filename ''.tsv''])');
  delete([filename '.tsv']);
  assert(isequal(c, [1;0;2]) && isequal(v, [1;1]) && isempty(dmu));
  assert(isequal(network3.actions, network.actions));
end
