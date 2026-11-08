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

%% kinetic regulation: network -> SBtab (column KineticRegulation) -> network
if has_mnt
  network.regulation_matrix = sparse([0 0 1; -1 0 0]);
  [~, d] = evalc('network_to_sbtab(network)');
  assert(isequal(sbtab_table_get_column(d.tables.Reaction,'KineticRegulation'), {'+ C'; '- A'}));
  assert(~sbtab_table_has_column(d.tables.Reaction,'MetabolicRegulation'));
  [~, network4] = evalc('sbtab_to_network(d)');
  assert(isequal(full(network4.regulation_matrix), [0 0 1; -1 0 0]));
end

%% kinetic regulation: older column names and notations
if has_mnt
  columns = struct('ID', {{'R1';'R2'}}, 'ReactionFormula', {{'A <=> B';'B <=> C'}});
  for my_case = {{'KineticRegulation', {'+ C'; '- A'}}, {'MetabolicRegulation', {'+ C'; '- A'}}, ...
                 {'Regulator', {'+C'; '-A'}}}
    my_columns = columns;
    my_columns.(my_case{1}{1}) = my_case{1}{2};
    [~, network5] = evalc('network_build_from_sum_formulae([], [], my_columns)');
    assert(isequal(full(network5.regulation_matrix), [0 0 1; -1 0 0]));
  end
  % regulators with dual or unspecified effect are ignored (with a warning)
  my_columns = columns;
  my_columns.Regulator = {'+C|B'; '+- A'};
  [~, network6] = evalc('network_build_from_sum_formulae([], [], my_columns)');
  assert(isequal(full(network6.regulation_matrix), [0 0 1; 0 0 0]));
end

%% kinetic regulation: unknown compounds give an error
if has_mnt
  my_columns = struct('ID', {{'R1';'R2'}}, 'ReactionFormula', {{'A <=> B';'B <=> C'}}, ...
                      'KineticRegulation', {{'+ X'; ''}});
  try
    evalc('network_build_from_sum_formulae([], [], my_columns)');
    assert(false, 'no error for an unknown regulator');
  catch err
    assert(contains(err.message, 'Unknown compound(s) X'));
  end
end

%% kinetic regulation: Regulation table (TableID 'Regulation') added to column KineticRegulation
if has_mnt
  doc_text = @(reaction_extra, regulation) sprintf(['!!!SBtab Document=''reg'' SBtabVersion=''1.2''\n' ...
      '!!SBtab TableID=''Compound'' TableType=''Compound''\n!ID\nA\nB\nC\n\n' ...
      '!!SBtab TableID=''Reaction'' TableType=''Reaction''\n!ID\t!ReactionFormula%s\n%s\n\n%s'], ...
      reaction_extra{:}, regulation);
  cases = { ...
    % new columns; RegulationType, RegulationType:SBO and RegulatoryCoefficient
    {{sprintf('\t!KineticRegulation'), sprintf('R1\tA <=> B\t+ C\nR2\tB <=> C\t')}, ...
     sprintf(['!!SBtab TableID=''Regulation'' TableType=''Regulation''\n' ...
              '!ID\t!Regulator\t!Target\t!RegulationType\t!RegulationType:SBO\t!RegulatoryCoefficient\n' ...
              'Reg1\tA\tR2\tinhibition\tSBO:0000020\t-2\n']), [0 0 1; -1 0 0], 1}, ...
    % columns RegulatorType and TargetType (other types are skipped); NaN: unknown or dual
    {{'', sprintf('R1\tA <=> B\nR2\tB <=> C')}, ...
     sprintf(['!!SBtab TableID=''Regulation'' TableType=''Regulation''\n' ...
              '!ID\t!Regulator\t!Target\t!RegulatorType\t!TargetType\t!RegulatoryCoefficient\n' ...
              'Reg1\tC\tR1\tmetabolite\tenzyme\t1\nReg2\tA\tR2\tmetabolite\tenzyme\tNaN\n' ...
              'Reg3\tTF1\tR2\ttranscriptional regulator\tpromoter\t-1\n']), [0 0 1; 0 0 0], 0}, ...
    % sign from the SBO term or from the coefficient only
    {{'', sprintf('R1\tA <=> B\nR2\tB <=> C')}, ...
     sprintf(['!!SBtab TableID=''Regulation'' TableType=''Regulation''\n' ...
              '!TargetReaction\t!Regulator\t!RegulationType:SBO\t!RegulatoryCoefficient\n' ...
              'R1\tC\tSBO:0000459\t\nR2\tA\t\t-1\n']), [0 0 1; -1 0 0], 0}, ...
    % older files: TableID 'EnzymeRegulation', columns Reaction and Compound
    {{'', sprintf('R1\tA <=> B\nR2\tB <=> C')}, ...
     sprintf(['!!SBtab TableID=''EnzymeRegulation'' TableType=''Regulation''\n' ...
              '!RegulationType\t!Reaction\t!Compound\nactivation\tR1\tC\ninhibition\tR2\tA\n']), [0 0 1; -1 0 0], 0}};
  for it = 1:length(cases)
    my_case = cases{it};
    filename = [tempname '.tsv'];
    fid = fopen(filename,'w'); fprintf(fid, '%s', doc_text(my_case{1}, my_case{2})); fclose(fid);
    d = sbtab_document_load_from_one(filename);
    delete(filename);
    lastwarn('');
    [~, network7] = evalc('sbtab_to_network(d, struct(''load_quantity_table'',0))');
    assert(isequal(full(network7.regulation_matrix), my_case{3}), sprintf('case %d', it));
    if my_case{4}
      assert(contains(lastwarn, 'both in the Reaction column'), sprintf('case %d: no warning', it));
    end
  end
end

%% regulation formulae (columns KineticRegulation and ActivityCondition)
if has_mnt
  [s, c] = network_parse_regulation_formula('+ ADP - PEP +- ATP');
  assert(isequal(s, [1 -1 0]) && isequal(c, {'ADP';'PEP';'ATP'}));
  [s, c] = network_parse_regulation_formula('-PEP|+ADP');
  assert(isequal(s, [-1 1]) && isequal(c, {'PEP';'ADP'}));
  [s, c] = network_parse_regulation_formula('');
  assert(isempty(s) && isempty(c));
end

%% Regulator table with column ActivityCondition: SBtab -> network -> SBtab
if has_mnt
  text = sprintf(['!!!SBtab Document=''reg'' SBtabVersion=''1.1''\n' ...
      '!!SBtab TableID=''Compound'' TableType=''Compound''\n!ID\nA\nB\nC\n\n' ...
      '!!SBtab TableID=''Reaction'' TableType=''Reaction''\n!ID\t!ReactionFormula\nR1\tA <=> B\nR2\tB <=> C\n\n' ...
      '!!SBtab TableID=''Regulator'' TableType=''Regulator''\n!ID\t!ActivityCondition\nTF1\t+ A\nTF2\t- C + B\nTF3\t\n']);
  filename = [tempname '.tsv'];
  fid = fopen(filename,'w'); fprintf(fid, '%s', text); fclose(fid);
  d = sbtab_document_load_from_one(filename);
  delete(filename);
  [~, network8] = evalc('sbtab_to_network(d, struct(''load_quantity_table'',0))');
  assert(isequal(network8.regulators, {'TF1';'TF2';'TF3'}));
  assert(isequal(full(network8.regulator_activity_matrix), [1 0 0; 0 1 -1; 0 0 0]));
  [~, d2] = evalc('network_to_sbtab(network8)');
  assert(isfield(d2.tables, 'Regulator'));
  assert(isequal(sbtab_table_get_column(d2.tables.Regulator,'ActivityCondition'), {'+ A';'+ B - C';''}));
  % unknown compounds are ignored with a warning
  network9 = sbtab_regulator_table_to_network(network8, sbtab_table_construct(struct('TableID','Regulator','TableType','Regulator'), {'ID','ActivityCondition'}, {{'TF1'},{'+ X'}}));
  assert(nnz(network9.regulator_activity_matrix) == 0);
  assert(contains(lastwarn, 'unknown compound'));
end
