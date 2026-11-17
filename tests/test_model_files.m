% Tests for model files: compartments in SBtab, mnt_save_model_and_state / mnt_load_model_and_state,
% layout table (netgraph_layout_table, netgraph_layout_from_table), network_choose, cobra_to_network

setup_test_path;

network = example_network;
network.metabolites = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}';
network.compartments            = {'cell'; 'medium'};
network.compartment_sizes       = [1; 2];
network.metabolite_compartments = {'medium', 'medium', 'cell', 'cell', 'medium', 'medium'}';
c   = [2 1 0.5 0.4 0.1 0.2]';
v   = [1 1 2 1.5 0.5]';
keq = [10 20 5 3 4]';

%% compartments: SBtab round trip (Compound!Location, table Compartment)
filename = [tempname '.tsv'];
sbtab_network_state_save(network, c, v, filename, [], struct('verbose', 0), keq);
network2 = sbtab_network_state_load(filename);
delete(filename);
assert(isequal(network2.metabolite_compartments(:), network.metabolite_compartments(:)));
assert(isequal(network2.compartments(:), network.compartments(:)));
assert(isequal(network2.compartment_sizes(:), network.compartment_sizes(:)));

%% network_choose keeps compartments and restricts biomass reactions
network.biomass_reactions = {'E5'; 'E3'};
sub = network_choose(network, 1:6, [1 2 3 4], 1);
assert(isequal(sub.compartments, network.compartments) && isequal(sub.compartment_sizes, network.compartment_sizes));
assert(isequal(sub.biomass_reactions, {'E3'}));
sub = network_remove_action(network, 5);
assert(isequal(sub.compartments, network.compartments));

%% mnt_save_model_and_state / mnt_load_model_and_state with layout table
networkG = netgraph_make_layout_network(network, {'S6'}, sbtab_table_construct(struct('TableID','Position','TableType','Position'), ...
             {'Element','Position:x','Position:y'}, {{'S1';'S2';'S3';'S4';'S5';'E1';'E2';'E3';'E4';'E5'}, (1:10)', (11:20)'}));
assert(length(networkG.metabolites) == 5);
outdir = tempname; mkdir(outdir);
mnt_save_model_and_state(fullfile(outdir, 'Model'), network, v, c, keq, networkG, struct('formats', {{'mat','tsv'}}));
[n_tsv, v_tsv, c_tsv, keq_tsv, ~, G_tsv] = mnt_load_model_and_state(fullfile(outdir, 'Model.tsv'));
[n_mat, v_mat, ~, ~, ~, G_mat] = mnt_load_model_and_state(fullfile(outdir, 'Model.mat'));
assert(max(abs(v_tsv - v)) < 1e-12 && max(abs(c_tsv - c)) < 1e-12 && max(abs(keq_tsv - keq)) < 1e-10);
assert(isequal(G_tsv.metabolites(:), networkG.metabolites(:)));
assert(max(max(abs(G_tsv.graphics_par.x - networkG.graphics_par.x))) < 1e-12);
assert(isequal(n_tsv.metabolite_compartments(:), network.metabolite_compartments(:)));
assert(isequal(G_mat.metabolites, networkG.metabolites) && isequal(v_mat, v));
% network without state
mnt_save_model_and_state(fullfile(outdir, 'NetworkModel'), network, [], [], [], networkG, struct('formats', {{'tsv'}}));
[n0, v0, c0, keq0, theta0, G0] = mnt_load_model_and_state(fullfile(outdir, 'NetworkModel.tsv'));
assert(isequal(full(n0.N), full(network.N)) && isempty(keq0) && isempty(theta0));
assert(~isempty(G0));
rmdir(outdir, 's');

%% cobra_to_network: exchange reactions by prefix, external metabolites by compartment, biomass from objective
model.mets     = {'a[e]'; 'a[c]'; 'b-L[c]'; 'c[c]'};
model.rxns     = {'EX_a'; 'Ta'; 'R1'; 'BIOMASS'};
model.S        = sparse([-1 -1 0 0; 0 1 -1 0; 0 0 1 -1; 0 0 0 1]);
model.rev      = [1; 1; 0; 0];
model.c        = [0; 0; 0; 1];
model.metNames = model.mets;
model.rxnNames = model.rxns;
n = cobra_to_network(model, struct('remove_reactions_prefix', 'EX_', 'external_method', 'compartment', ...
      'external_single_reaction', 1, 'metabolite_id_replacements', {{'-L', ''}}, 'bigg_ids', 1, 'full_matrix', 1));
assert(isequal(n.actions, {'Ta'; 'R1'; 'BIOMASS'}));
assert(isequal(n.metabolites, {'a[e]'; 'a[c]'; 'b[c]'; 'c[c]'}));
assert(isequal(n.external(:)', [1 0 0 1]));          % a[e]: compartment e; c[c]: only in one reaction
assert(isequal(n.biomass_reactions, {'BIOMASS'}));
assert(isequal(n.metabolite_BIGGID, n.metabolites));
n_default = cobra_to_network(model);                  % default: as before 2026-11-17
assert(length(n_default.actions) == 4 && isequal(n_default.reversible(:)', [1 1 0 0]));

%% network_add_reaction keeps all fields
network = example_network;
network.metabolite_names        = network.metabolites;
network.metabolite_compartments = repmat({'cell'}, 6, 1);
network.compartments = {'cell'}; network.compartment_sizes = 1;
network.reaction_names  = network.actions;
network.reaction_KEGGID = repmat({'R00001'}, 5, 1);
network.enzyme_mass     = ones(5,1);
n2 = network_add_reaction(network, {'S3','S4'}, [-1 1], {'S3/S4 conversion'}, 1);
assert(isequal(n2.actions{end}, 'S3/S4 conversion') && isequal(n2.reaction_names{end}, 'S3/S4 conversion'));
assert(isequal(n2.reaction_KEGGID{end}, '') && isnan(n2.enzyme_mass(end)) && length(n2.reversible) == 6);
assert(isequal(n2.metabolite_compartments, network.metabolite_compartments) && isequal(n2.compartments, {'cell'}));
assert(isequal(full(n2.N(:,end))', [0 0 -1 1 0 0]) && isequal(size(n2.regulation_matrix), [6 6]));
n3 = network_add_reaction(network, {'S3','X'}, [-1 1], 'new');   % new metabolite X (external)
assert(length(n3.metabolites) == 7 && n3.external(7) == 1 && length(n3.metabolite_compartments) == 7);
