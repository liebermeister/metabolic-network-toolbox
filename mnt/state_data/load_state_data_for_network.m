function [M, ids, states] = load_state_data_for_network(network, data_file, quantity_type, file_options, table_id)

% [M, ids, states] = load_state_data_for_network(network, data_file, quantity_type, file_options)
%
% Load flux, metabolite, or enzyme data from SBtab data file
%
% quantity_type: 'metabolite_concentration', 'enzyme_concentration', 'reaction_rate', 'reaction_gibbs_free_energy'
%
% default file_options values:
%   file_options.columns_mean           = [];
%   file_options.columns_std            = [];
%   file_options.replace_ids_in_network = [];
%   file_options.match_data_by          = 'ModelElementId'; % KeggId; 'BiggId'
%   file_options.verbose                = 0;
%
% Metabolite data matched by 'BiggId': the BiGG ids of the data (column Compound:Identifiers:bigg.compound)
% are matched to network.metabolite_BIGGID by bigg_match_metabolite_ids, so compartments may be written
% differently ('glc__D', 'glc__D_c', 'glc__D[c]'; '' and 'c' count as the same compartment). The mapping
% is printed if file_options.verbose = 1 or if there are problems (data ids not found, ambiguous, or
% mapped with a different compartment). (Until 2026-11-15, BiGG ids had to be identical.)
% In all other cases, the data ids must be identical to the ids in the network.
%
% To load a model (with kinetic data) and state data, use sbtab_load_network_model_and_data

if isempty(data_file)
  M = []; ids= []; states = [];
  return
end

eval(default('file_options','struct', 'table_id','[]'));
file_options_default.columns_mean = [];
file_options_default.columns_std = [];
file_options_default.replace_ids_in_network = [];
file_options_default.match_data_by = 'ModelElementId'; % 'KeggId'; ; 'BiggId'%
file_options_default.verbose = 0;
file_options = join_struct(file_options_default,file_options);

if isempty(file_options.columns_mean)
  error('No column names specified');
end

mean_columns = column(file_options.columns_mean);
std_columns  = column(file_options.columns_std);
replace_ids_in_network = file_options.replace_ids_in_network;

switch quantity_type

  case 'metabolite_concentration'
    quantity_type = 'concentration';
    table_name = 'MetaboliteConcentration';
    switch file_options.match_data_by
      case 'BiggId'
        element_IDs   = network.metabolite_BIGGID;
        id_column     = '!Compound:Identifiers:bigg.compound';
      case 'KeggId'
        element_IDs   = network.metabolite_KEGGID;
        id_column     = '!Compound:Identifiers:kegg.compound';
      case 'ModelElementId'
        element_IDs   = network.metabolites;
        id_column     = '!Compound';
    end

  case 'enzyme_concentration'
    quantity_type = 'concentration of enzyme';
    table_name = 'EnzymeConcentration';
    switch file_options.match_data_by
      case 'BiggId'
        element_IDs   = network.reaction_BIGGID;
        id_column     = '!Reaction:Identifiers:bigg.reaction';
      case 'KeggId'
        element_IDs   = network.reaction_KEGGID;
        id_column     = '!Reaction:Identifiers:kegg.reaction';
      case 'ModelElementId'
        element_IDs   = network.actions;
        id_column     = '!Reaction';
    end

  case 'reaction_rate'
    quantity_type = 'rate of reaction';
    table_name = 'Flux';
    switch file_options.match_data_by
      case 'BiggId'
        element_IDs   = network.reaction_BIGGID;
        id_column     = '!Reaction:Identifiers:bigg.reaction';
      case 'KeggId'
        element_IDs   = network.reaction_KEGGID;
        id_column     = '!Reaction:Identifiers:kegg.reaction';
      case 'ModelElementId'
        element_IDs   = network.actions;
        id_column     = '!Reaction';
    end

  case 'reaction_gibbs_free_energy'
    quantity_type = 'Gibbs energy of reaction';
    table_name = 'ReactionGibbsFreeEnergy';
    switch file_options.match_data_by
      case 'BiggId'
        element_IDs   = network.reaction_BIGGID;
        id_column     = '!Reaction:Identifiers:bigg.reaction';
      case 'KeggId'
        element_IDs   = network.reaction_KEGGID;
        id_column     = '!Reaction:Identifiers:kegg.reaction';
      case 'ModelElementId'
        element_IDs   = network.actions;
        id_column     = '!Reaction';
    end

  otherwise
    error(sprintf('unknown data type %s',quantity_type));

end

if length(replace_ids_in_network)
  for it = 1:length(replace_ids_in_network)
    replace = replace_ids_in_network{it}{1};
    by      = replace_ids_in_network{it}{2};
    ind = label_names(replace,element_IDs);
    if ind
      element_IDs{ind} = by;
    end
  end
end

if length(table_id), table_name = table_id; end

% metabolite BiGG ids: match data ids to model ids (compartment notations may differ)
model_IDs = element_IDs;
if strcmp(quantity_type, 'concentration') && strcmp(file_options.match_data_by, 'BiggId')
  if ischar(data_file), data_file = sbtab_document_load(data_file); end
  [~, data_IDs] = sbtab_load_quantity_data(data_file, table_name, quantity_type, id_column, [], mean_columns(1), 0);
  [ind_match, match_report] = bigg_match_metabolite_ids(element_IDs, data_IDs, network.external);
  element_IDs = repmat({''}, size(element_IDs));
  element_IDs(ind_match > 0) = data_IDs(ind_match(ind_match > 0));
  if file_options.verbose || match_report.has_problems
    bigg_match_metabolite_ids_display(match_report, 'Metabolite data (BiGG ids) mapped to the model');
  end
end

[M.Mean, ~, states, M.Unit] = sbtab_load_quantity_data(data_file, table_name, quantity_type, id_column, element_IDs, mean_columns,1);
ids = model_IDs;

if length(std_columns)
  M.Std = sbtab_load_quantity_data(data_file, [], quantity_type, id_column, element_IDs, std_columns,1);
else
  M.Std = [];
end

M.ids = ids;

% units: concentrations (metabolites, enzymes) are converted to mM and Gibbs free energies to kJ/mol,
% the units used in MNT; unknown units give an error
if ~isempty(M.Unit)
  f = 1; unit_out = M.Unit;
  switch quantity_type
    case {'concentration', 'concentration of enzyme'}
      f = concentration_unit_in_mM(M.Unit, sprintf('quantity "%s", table %s', quantity_type, table_name));
      unit_out = 'mM';
    case 'Gibbs energy of reaction'
      switch strrep(M.Unit, ' ', '')
        case 'kJ/mol', f = 1;
        case 'J/mol',  f = 0.001;
        otherwise, error('Unknown unit ''%s'' for quantity "%s" (allowed: kJ/mol, J/mol)', M.Unit, quantity_type);
      end
      unit_out = 'kJ/mol';
  end
  if f ~= 1
    M.Mean = f * M.Mean;
    if ~isempty(M.Std), M.Std = f * M.Std; end
    display(sprintf('- (load_state_data_for_network.m): converting %s data from %s to %s', quantity_type, M.Unit, unit_out));
  end
  M.Unit = unit_out;
end

if prod(isnan(M.Mean))
  warning(sprintf('No data for data type "%s" found in data file',quantity_type))
end
