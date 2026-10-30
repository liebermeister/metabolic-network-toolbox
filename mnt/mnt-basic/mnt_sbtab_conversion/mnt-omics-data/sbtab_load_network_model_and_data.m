function result = sbtab_load_network_model_and_data(data_file, file_options)

% result = sbtab_load_network_model_and_data(data_file, file_options)
%
% Inputs
%   data_file: (string or cell array of strings) SBtab file(s) containing tables with input data,
%              describing a metabolic model, kinetic model parameters (in table "Parameter"), and metabolic state data
%
%   file_options:   struct with (optional) fields describing details of the state data file
%      file_options.match_data_by       'ModelElementId' or 'KeggId' or 'BiggId';
%      file_options.flux_table_id       (default 'FluxData')
%      file_options.metabolite_table_id (default 'MetaboliteConcentrationData')
%      file_options.enzyme_table_id     (default 'EnzymeConcentrationData')
%      file_options.delta_g_table_id    (default 'ReactionGibbsFreeEnergyData')
%      file_options.columns_mean:       cell array of column names for mean data values
%                                       (same columns names in metabolite, flux, and enzyme table!)
%      file_options.columns_std:        same, for std dev columns
%   also:
%      file_options.network             network model (needed if the input files do not contain the network model)
%      file_options.load_quantity_table generate field "network.kinetics" from Parameter table (if it exists) (default 0)
%      file_options.prior_file          prior file for parameter_table (optional, default value = []; then cmb_prior_file is used)
%
% This function requires the packages Parameter balancing (function parameter_balancing_prior)
% and, if no prior file is given, Model balancing (function cmb_prior_file)
%
% Output
%   result.network
%   result.kinetic_data
%   result.state_data
%   result.file_options
%
% For saving network model (including kinetic data) and state data, see: sbtab_save_network_model_and_data

eval(default('file_options','struct'));

file_options_default.metabolite_table_id = 'MetaboliteConcentrationData';
file_options_default.flux_table_id       = 'FluxData';
file_options_default.enzyme_table_id     = 'EnzymeConcentrationData';
file_options_default.delta_g_table_id    = 'ReactionGibbsFreeEnergyData';
file_options_default.columns_mean        = {'Mean'};
file_options_default.columns_std         = {'Std'};
file_options_default.match_data_by       = 'ModelElementId'; % 'KeggId'; % 'BiggId'
file_options_default.load_quantity_table = 0;
file_options_default.prior_file          = [];
file_options = join_struct(file_options_default, file_options);

% load SBtab document from one or several files

if iscell(data_file)
  s = sbtab_document_load(data_file);
else
  s = sbtab_document_load_from_one(data_file);
end

% extract network

if label_names({'Reaction'},sbtab_document_get_table_names(s))
  opt = struct('load_quantity_table',file_options.load_quantity_table);
  result.network = sbtab_to_network(s, opt);
else
  display('No network information found in file');
  if isfield(file_options,'network')
    result.network = file_options.network;
  else
    error('Please provide network structure in file_options.network');
  end
end

% extract kinetic constants (if available)
if isempty(file_options.prior_file)
  file_options.prior_file = cmb_prior_file;
end

parameter_prior = parameter_balancing_prior([], file_options.prior_file, 0);

if label_names('ParameterData', sbtab_document_get_table_names(s))
  result.kinetic_data = kinetic_data_load([],parameter_prior, result.network, sbtab_document_get_table(s,'ParameterData'), struct('enforce_ranges',0));
elseif label_names('Parameter', sbtab_document_get_table_names(s))
  result.kinetic_data = kinetic_data_load([],parameter_prior, result.network, sbtab_document_get_table(s,'Parameter'), struct('enforce_ranges',0));
else
  warning('No kinetic data table (ID "Parameter") found');
end

% check
% kinetic_data_print(result.kinetic_data,result.network)

% extract metabolic state data

if length(file_options.metabolite_table_id)
  data_metabolite = sbtab_document_get_table(s,file_options.metabolite_table_id);
  [result.state_data.metabolite_data, ids, states] = load_state_data_for_network(result.network, data_metabolite, 'metabolite_concentration', file_options);
end

if length(file_options.flux_table_id)
  data_flux = sbtab_document_get_table(s,file_options.flux_table_id);
  [result.state_data.flux_data, ids, states] = load_state_data_for_network(result.network, data_flux, 'reaction_rate', file_options);
end

if length(file_options.enzyme_table_id)
  data_enzyme     = sbtab_document_get_table(s,file_options.enzyme_table_id);
  [result.state_data.enzyme_data, ids, states] = load_state_data_for_network(result.network, data_enzyme, 'enzyme_concentration', file_options);
end

if length(file_options.delta_g_table_id)
  data_delta_g    = sbtab_document_get_table(s,file_options.delta_g_table_id);
  [result.state_data.delta_g_data, ids, states] = load_state_data_for_network(result.network, data_delta_g, 'reaction_gibbs_free_energy', file_options);
end

% extract tool options (if available)

if label_names('Config', sbtab_document_get_table_names(s))
  result.file_options    = sbtab_to_options(sbtab_document_get_table(s,'Config'));
else
  %  warning('No options table (ID "Config") found');
end
