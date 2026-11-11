function sbml_object = sbtab_to_sbml(sbtab_object, filename, options)

% sbml_object = sbtab_to_sbml(sbtab_object, filename, options)
%
% Convert an SBtab model (SBtab document or file name) into an SBML model
% (Matlab struct in libSBML's format). If 'filename' is given, the model is
% also written to an SBML file.
%
% Metabolic state: if the SBtab document contains state tables (TableIDs MetaboliteConcentration,
% Flux, ReactionGibbsFreeEnergy, EquilibriumConstant), the model is written with
% sbml_network_state_save, in the SBML format of the Metabolic Network Toolbox for models with a
% metabolic state: concentrations as initial concentrations of the species; fluxes, equilibrium
% constants and reaction Gibbs free energies as local parameters 'reference_flux',
% 'equilibrium_constant', 'reaction_gibbs_free_energy' of the reactions (with units, SBO terms
% and notes; see sbml_network_state_save). Otherwise, the network is written with network_sbml_export.
%
% options (optional struct)
%   .state   'auto' (default: write the state if the document contains state tables),
%            1 (always write the state format), 0 (never)
%   other fields are passed to sbml_network_state_save (e.g. .name, .concentration_unit, .time_unit)
%
% This function requires the Metabolic Network Toolbox and libSBML (see network_sbml_export)
%
% usage example:
% sbtab_object = sbtab_load_example('sbtab/escherichia-coli/escherichia-coli_noor_2016_model.tsv','sbtab');
% sbml_object = sbtab_to_sbml(sbtab_object);
%
% See also: sbml_to_sbtab, sbml_network_state_save, sbtab_network_state_load

eval(default('filename','[]','options','struct'));
if isempty(options), options = struct; end
options = join_struct(struct('state','auto','verbose',0), options);

sbtab_object = sbtab_get(sbtab_object);

state_tables = {'MetaboliteConcentration', 'Flux', 'ReactionGibbsFreeEnergy', 'EquilibriumConstant'};
switch options.state
  case 'auto', write_state = any(ismember(state_tables, sbtab_document_get_table_names(sbtab_object)));
  otherwise,   write_state = options.state;
end

if write_state
  warning_state = warning('off', 'all');   % missing state tables are fine here
  [network, c, v, dmu, keq] = sbtab_network_state_load(sbtab_object);
  warning(warning_state);
  save_options = rmfield(options, 'state');
  sbml_object = sbml_network_state_save(network, c, v, filename, dmu, save_options, keq);
else
  sbml_object = network_sbml_export(sbtab_to_network(sbtab_object), 0, 'Model', filename);
end
