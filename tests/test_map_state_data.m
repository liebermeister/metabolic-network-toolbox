function tests = test_map_state_data

% Tests for mnt_map_state_data (mapping of flux and concentration data from SBtab files to a model;
% needs the SBtab toolbox)
%
% Run:  results = runtests('test_map_state_data')

tests = functiontests(localfunctions);


% ------------------------------------------------------------------
% small test network: R1: A <=> B, R2: B + C <=> D, R3: D <=> E

function network = make_network()

network.actions     = {'R1'; 'R2'; 'R3'};
network.metabolites = {'A'; 'B'; 'C'; 'D'; 'E'};
network.metabolite_KEGGID = {'C1'; 'C2'; 'C3'; 'C4'; 'C5'};
network.reaction_KEGGID   = {'K1'; 'K2'; 'K3'};
network.N = [-1 0 0; 1 -1 0; 0 -1 0; 0 1 -1; 0 0 1];
network.external = [1; 0; 1; 0; 1];


function filename = write_data_file(testCase, flux_formula_R2)

% data file with two conditions (columns x_Mean, y_Mean); no data for R3 and for metabolite E;
% the reaction formula of R2 is given in the data file (to test the orientation)
filename = [tempname '.tsv'];
fid = fopen(filename, 'w');
fprintf(fid, '!!!SBtab Document=''test''\n\n');
fprintf(fid, '!!SBtab TableID=''FluxData'' TableType=''QuantityMatrix'' Unit=''mM/s''\n');
fprintf(fid, '!QuantityType\t!Reaction\t!ReactionFormula\t!Reaction:Identifiers:kegg.reaction\tx_Mean\tx_Std\ty_Mean\ty_Std\n');
fprintf(fid, 'rate of reaction\tR1\tA <=> B\tK1\t1\t0.1\t2\t0.2\n');
fprintf(fid, 'rate of reaction\tR2\t%s\tK2\t0.5\t0.05\t-1\t0.1\n', flux_formula_R2);
fprintf(fid, '\n!!SBtab TableID=''MetaboliteConcentrationData'' TableType=''QuantityMatrix'' Unit=''mM''\n');
fprintf(fid, '!QuantityType\t!Compound\t!Compound:Identifiers:kegg.compound\tx_Mean\tx_Std\ty_Mean\ty_Std\n');
fprintf(fid, 'concentration\tA\tC1\t1\t0.1\t1.5\t0.1\n');
fprintf(fid, 'concentration\tB\tC2\t2\t0.2\t2.5\t0.2\n');
fprintf(fid, 'concentration\tC\tC3\t3\t0.3\t3.5\t0.3\n');
fprintf(fid, 'concentration\tD\tC4\t4\t0.4\t4.5\t0.4\n');
fclose(fid);
testCase.addTeardown(@() delete(filename));


% ------------------------------------------------------------------

function test_mapping_by_model_element_ids(testCase)

network = make_network;
file = write_data_file(testCase, 'B + C <=> D');
data = mnt_map_state_data(network, file, struct('columns_mean', {{'x_Mean', 'y_Mean'}}, 'columns_std', {{'x_Std', 'y_Std'}}));
verifyEqual(testCase, data.v_data, [1 2; 0.5 -1; nan nan]);
verifyEqual(testCase, data.v_std, [0.1 0.2; 0.05 0.1; nan nan]);
verifyEqual(testCase, data.c_data, [1 1.5; 2 2.5; 3 3.5; 4 4.5; nan nan]);
verifyEmpty(testCase, data.e_data);


function test_mapping_by_kegg_ids(testCase)

network = make_network;
file = write_data_file(testCase, 'B + C <=> D');
data = mnt_map_state_data(network, file, struct('match_data_by', 'KeggId', 'columns_mean', {{'y_Mean'}}));
verifyEqual(testCase, data.v_data, [2; -1; nan]);
verifyEqual(testCase, data.c_data, [1.5; 2.5; 3.5; 4.5; nan]);


function test_reaction_orientation(testCase)

% data for R2 given in the opposite orientation: the sign is changed (orient_data = 1, default),
% or kept (orient_data = 0)
network = make_network;
file = write_data_file(testCase, 'D <=> B + C');
data = mnt_map_state_data(network, file, struct('columns_mean', {{'x_Mean'}}));
verifyEqual(testCase, data.v_data, [1; -0.5; nan]);
data = mnt_map_state_data(network, file, struct('columns_mean', {{'x_Mean'}}, 'orient_data', 0));
verifyEqual(testCase, data.v_data, [1; 0.5; nan]);


function test_selected_tables(testCase)

network = make_network;
file = write_data_file(testCase, 'B + C <=> D');
data = mnt_map_state_data(network, file, struct('columns_mean', {{'x_Mean'}}, 'metabolite_table_id', []));
verifyEmpty(testCase, data.c_data);
verifyEqual(testCase, data.v_data, [1; 0.5; nan]);
verifyError(testCase, @() mnt_map_state_data(network, file, struct('enzyme_table_id', 'EnzymeConcentrationData')), ?MException);
