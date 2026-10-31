function sbtab_model_display(sbtab_file)

% sbtab_model_display(sbtab_file)
%
% Read SBtab toolkit SBtab file and print / plot some overview information
%
% Usage example:
% dir = [sbtab_basedir '/../resources/examples/sbtab/ecoli/'];
% sbtab_model_display([dir 'ecoli_noor_2016_model.tsv']);

sbtab_object = sbtab_document_load(sbtab_file);

sbtab_document_print(sbtab_object)

[network, kinetic_data] = sbtab_to_network(sbtab_file);

netgraph_concentrations(network,[],[],1);

%kinetic_data_print(kinetic_data, network)

kinetics_print(network.kinetics,1)

%  [v, v_std] = flux_data_load(network, flux_data_file, reaction_field, reaction_column, value_column, std_column)
