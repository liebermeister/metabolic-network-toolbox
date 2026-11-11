function tests = test_network_add_kinetic_regulation

% Tests for network_add_kinetic_regulation (kinetic regulation from an SBtab Reaction table)
%
% Run:  results = runtests('test_network_add_kinetic_regulation')

tests = functiontests(localfunctions);


function network = make_network()

network.actions     = {'R1'; 'R2'};
network.metabolites = {'A'; 'B'; 'C'};
network.N = [-1 0; 1 -1; 0 1];
network.external = [1; 0; 1];
network.regulation_matrix = zeros(2, 3);


function filename = write_file(testCase, rows)

filename = [tempname '.tsv'];
fid = fopen(filename, 'w');
fprintf(fid, '!!!SBtab Document=''test''\n%% comment line\n\n!!SBtab TableID=''Reaction'' TableType=''Reaction''\n!ID\t!KineticRegulation\n');
for it = 1:size(rows, 1), fprintf(fid, '%s\t%s\n', rows{it, 1}, rows{it, 2}); end
fclose(fid);
testCase.addTeardown(@() delete(filename));


function test_add_from_file(testCase)

network = make_network;
network.regulation_matrix(1, 3) = 1;
network = network_add_kinetic_regulation(network, write_file(testCase, {'R2', '- A'; 'R1', '+ B - C'}));
verifyEqual(testCase, full(network.regulation_matrix), [0 1 -1; -1 0 0]);


function test_replace(testCase)

network = make_network;
network.regulation_matrix(1, 1) = 1;
network.regulation_matrix(2, 1) = 1;
network = network_add_kinetic_regulation(network, write_file(testCase, {'R1', '- C'}), struct('replace', 1));
verifyEqual(testCase, full(network.regulation_matrix), [0 0 -1; 1 0 0]);


function test_errors(testCase)

network = make_network;
verifyError(testCase, @() network_add_kinetic_regulation(network, write_file(testCase, {'R9', '- A'})), ?MException);
verifyError(testCase, @() network_add_kinetic_regulation(network, write_file(testCase, {'R1', '- X'})), ?MException);
