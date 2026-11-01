function sbml_object = sbtab_to_sbml(sbtab_object, filename)

% sbml_object = sbtab_to_sbml(sbtab_object, filename)
%
% Convert an SBtab model (SBtab document or file name) into an SBML model
% (Matlab struct in libSBML's format). If 'filename' is given, the model is
% also written to an SBML file.
%
% This function requires the Metabolic Network Toolbox and libSBML (see network_sbml_export)
%
% usage example:
% sbtab_object = sbtab_load_example('sbtab/ecoli/ecoli_noor_2016_model.tsv','sbtab');
% sbml_object = sbtab_to_sbml(sbtab_object);

eval(default('filename','[]'));

sbml_object = network_sbml_export(sbtab_to_network(sbtab_get(sbtab_object)), 0, 'Model', filename);
