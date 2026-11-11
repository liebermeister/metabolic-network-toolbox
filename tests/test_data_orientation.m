function tests = test_data_orientation

% Tests for the conversion of data tables to the reaction orientation of a model
% (functions in mnt/data_orientation; needs the SBtab toolbox)
%
% Run:  results = runtests('test_data_orientation')

tests = functiontests(localfunctions);


% ------------------------------------------------------------------
% small test network: R1: A <=> B, R2: B + C <=> D (KEGG IDs C1 ... C4)

function network = make_network()

network.actions     = {'R1'; 'R2'};
network.metabolites = {'A'; 'B'; 'C'; 'D'};
network.metabolite_KEGGID = {'C1'; 'C2'; 'C3'; 'C4'};
network.reaction_KEGGID   = {'K1'; 'K2'};
network.N = [-1 0; 1 -1; 0 -1; 0 1];


function t = make_table(quantity, reaction, formula, formula_kegg, value, attributes)

if ~exist('attributes', 'var'), attributes = struct; end
attributes.TableID   = 'Data';
attributes.TableType = 'Quantity';
t = sbtab_table_construct(attributes, ...
      {'QuantityType', 'Reaction', 'ReactionFormula', 'ReactionFormula:Identifiers:kegg.compound', 'Value', 'Min', 'Max'}, ...
      {quantity, reaction, formula, formula_kegg, value, repmat({'NaN'}, length(value), 1), repmat({'NaN'}, length(value), 1)});


% ------------------------------------------------------------------

function test_orientation_from_formulas(testCase)

network = make_network;
t = make_table({'rate of reaction'; 'rate of reaction'; 'rate of reaction'; 'rate of reaction'}, ...
               {'R1'; 'R1'; 'R2'; 'R2'}, ...
               {'A <=> B'; 'B <=> A'; ''; 'D <=> B + C'}, ...
               {''; ''; 'C2 + C3 <=> C4'; ''}, ...
               {'1'; '1'; '1'; '1'});
[orientation, ind] = sbtab_reaction_orientation(network, t);
verifyEqual(testCase, orientation, [1; -1; 1; -1]);
verifyEqual(testCase, ind, [1; 1; 2; 2]);


function test_orientation_from_table_attribute(testCase)

network = make_network;
t = make_table({'rate of reaction'}, {'R1'}, {''}, {''}, {'1'}, struct('ReactionOrientation', 'Model:test'));
verifyEqual(testCase, sbtab_reaction_orientation(network, t, struct('model_orientation', 'Model:test')), 1);
verifyTrue(testCase, isnan(sbtab_reaction_orientation(network, t, struct('model_orientation', 'Model:other'))));
% formulas have priority over the table attribute
t = make_table({'rate of reaction'}, {'R1'}, {'B <=> A'}, {''}, {'1'}, struct('ReactionOrientation', 'Model:test'));
verifyEqual(testCase, sbtab_reaction_orientation(network, t, struct('model_orientation', 'Model:test')), -1);


function test_matching_by_kegg_ids(testCase)

network = make_network;
t = make_table({'rate of reaction'}, {'X'}, {''}, {'C2 <=> C1'}, {'1'});
t = sbtab_table_set_column(t, 'Reaction:Identifiers:kegg.reaction', {'K1'});
[orientation, ind] = sbtab_reaction_orientation(network, t, struct('match_data_by', 'KeggId'));
verifyEqual(testCase, [orientation, ind], [-1, 1]);


function test_conversion_of_values(testCase)

network = make_network;
t = make_table({'rate of reaction'; 'equilibrium constant'; 'substrate catalytic rate constant'; 'Michaelis constant'; 'rate of reaction'}, ...
               {'R1'; 'R1'; 'R1'; 'R1'; 'R2'}, ...
               {'B <=> A'; 'B <=> A'; 'B <=> A'; 'B <=> A'; 'B + C <=> D'}, ...
               {''; ''; ''; ''; ''}, ...
               {'2'; '4'; '10'; '0.5'; '3'});
t = sbtab_table_set_column(t, 'Min', {'1'; '2'; 'NaN'; 'NaN'; 'NaN'});
t = sbtab_table_set_column(t, 'Max', {'3'; '8'; 'NaN'; 'NaN'; 'NaN'});
[t2, report] = sbtab_orient_data_table(network, t, struct('model_orientation', 'Model:test'));
verifyEqual(testCase, sbtab_table_get_column(t2, 'Value', 1), [-2; 0.25; 10; 0.5; 3], 'AbsTol', 1e-12);
verifyEqual(testCase, sbtab_table_get_column(t2, 'Min', 1), [-3; 0.125; nan; nan; nan], 'AbsTol', 1e-12);
verifyEqual(testCase, sbtab_table_get_column(t2, 'Max', 1), [-1; 0.5; nan; nan; nan], 'AbsTol', 1e-12);
q = sbtab_table_get_column(t2, 'QuantityType');
verifyEqual(testCase, q{3}, 'product catalytic rate constant');
f = sbtab_table_get_column(t2, 'ReactionFormula');
verifyEqual(testCase, f{1}, 'A <=> B');
verifyEqual(testCase, t2.attributes.ReactionOrientation, 'Model:test');
verifyEqual(testCase, numel(report.reversed_reactions), 3);   % the Michaelis constant does not change
verifyEmpty(testCase, report.unknown_orientation);


function test_unknown_orientation_is_reported(testCase)

network = make_network;
t = make_table({'rate of reaction'; 'concentration of enzyme'}, {'R1'; 'R1'}, {''; ''}, {''; ''}, {'2'; '1'});
[t2, report] = sbtab_orient_data_table(network, t, struct('model_orientation', 'Model:test'));
verifyEqual(testCase, report.unknown_orientation, {'R1'});   % only the flux row (enzyme levels do not depend on the orientation)
verifyEqual(testCase, sbtab_table_get_column(t2, 'Value', 1), [2; 1]);
verifyFalse(testCase, isfield(t2.attributes, 'ReactionOrientation'));


function test_formulas_of_network(testCase)

[f, fk] = network_reaction_formulas(make_network);
verifyEqual(testCase, f,  {'A <=> B'; 'B + C <=> D'});
verifyEqual(testCase, fk, {'C1 <=> C2'; 'C2 + C3 <=> C4'});
