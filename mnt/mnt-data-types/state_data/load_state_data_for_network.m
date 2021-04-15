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
%   file_options.match_data_by          = 'ModelElementId'; % KeggId
% 
% To load a model (with kinetic data) and state data, use sbtab_load_network_model_and_data

eval(default('file_options','struct', 'table_id','[]'));
file_options_default.columns_mean = [];
file_options_default.columns_std = [];
file_options_default.replace_ids_in_network = [];
file_options_default.match_data_by = 'ModelElementId'; % 'KeggId'; % 
file_options = join_struct(file_options_default,file_options);
  
if isempty(file_options.columns_mean),
  error('No column names specified');
end

mean_columns = column(file_options.columns_mean);
std_columns  = column(file_options.columns_std);
replace_ids_in_network = file_options.replace_ids_in_network;

switch quantity_type,
  
  case 'metabolite_concentration',
    quantity_type = 'concentration';
    table_name = 'MetaboliteConcentration';
    switch file_options.match_data_by,
      case 'KeggId',
        element_IDs   = network.metabolite_KEGGID;
        id_column     = '!Compound:Identifiers:kegg.compound';
      case 'ModelElementId'
        element_IDs   = network.metabolites;
        id_column     = '!Compound';
    end
  
  case 'enzyme_concentration',
    quantity_type = 'concentration of enzyme';      
    table_name = 'EnzymeConcentration';
    switch file_options.match_data_by,
      case 'KeggId',
        element_IDs   = network.reaction_KEGGID;
        id_column     = '!Reaction:Identifiers:kegg.reaction';
      case 'ModelElementId'
        element_IDs   = network.actions;
        id_column     = '!Reaction';
    end
  
  case 'reaction_rate',
    quantity_type = 'rate of reaction'; 
    table_name = 'MetabolicFlux';
    switch file_options.match_data_by,
      case 'KeggId',
        element_IDs   = network.reaction_KEGGID;
        id_column     = '!Reaction:Identifiers:kegg.reaction';
      case 'ModelElementId'
        element_IDs   = network.actions;
        id_column     = '!Reaction';
    end
  
  case 'reaction_gibbs_free_energy',
    quantity_type = 'Gibbs energy of reaction'; 
    table_name = 'ReactionGibbsFreeEnergy';
    switch file_options.match_data_by,
      case 'KeggId',
        element_IDs   = network.reaction_KEGGID;
        id_column     = '!Reaction:Identifiers:kegg.reaction';
      case 'ModelElementId'
        element_IDs   = network.actions;
        id_column     = '!Reaction';
    end

  otherwise 
    error(sprintf('unknown data type %s',quantity_type));
    
end

if length(replace_ids_in_network),
  for it = 1:length(replace_ids_in_network),
    replace = replace_ids_in_network{it}{1};
    by      = replace_ids_in_network{it}{2};
    ind = label_names(replace,element_IDs);
    if ind,
      element_IDs{ind} = by;
    end
  end
end

if length(table_id), table_name = table_id; end

[M.Mean, ids, states, M.Unit] = sbtab_load_quantity_data(data_file, table_name, quantity_type, id_column, element_IDs, mean_columns,1);

if length(std_columns),
  M.Std = sbtab_load_quantity_data(data_file, [], quantity_type, id_column, element_IDs, std_columns,1);
else
  M.Std = [];
end

M.ids = ids;

if prod(isnan(M.Mean)),
  warning(sprintf('No data for data type "%s" found in data file',quantity_type))
end
