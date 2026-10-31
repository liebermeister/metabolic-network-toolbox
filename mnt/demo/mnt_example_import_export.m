% MNT_EXAMPLE_IMPORT_EXPORT Save a network to an SBtab file and load it again
%
% This example requires SBtab for Matlab. For SBML files, which require libSBML, use
% network_sbml_export and network_sbml_import; for COBRA models, see 'help io'.
% For the network, see mnt_example_network

N = [ -1  0  0  0  0; ...
       0 -1  0  0  0; ...
       1  1 -1  0  0; ...
       0  0  1 -1 -1; ...
       0  0  0  1  0; ...
       0  0  0  0  1];
reversible   = [1 1 1 1 1]';
ind_external = [1 2 5 6]';
network      = network_construct(N, reversible, ind_external);

% Save the network to a file in the folder for temporary files

filename = fullfile(tempdir, 'mnt_example_model.tsv');
network_to_sbtab(network, struct('filename', filename));

% Load the network

network_2 = sbtab_to_network(filename, struct('one_sbtab_file', 1));
network_print_formulae(network_2)

delete(filename);
