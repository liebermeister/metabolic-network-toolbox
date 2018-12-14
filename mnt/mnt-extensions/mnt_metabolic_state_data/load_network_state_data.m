function [M, ids, states] = load_network_state_data(network,data_file,quantity_type,mean_columns,std_columns)

eval(default('std_columns','[]'));
  
switch quantity_type,
  
  case 'metabolite_concentration',
    quantity_type = 'concentration';
    element_IDs   = network.metabolite_KEGGID;
    id_column     = '!Compound:Identifiers:kegg.compound';
  
  case 'enzyme_concentration',
    quantity_type = 'concentration of enzyme';      
    element_IDs   = network.reaction_KEGGID;
    id_column     = '!Reaction:Identifiers:kegg.reaction';
  
  case 'reaction_flux',
    quantity_type = 'flux';      
    element_IDs   = network.reaction_KEGGID;
    id_column     = '!Reaction:Identifiers:kegg.reaction';
    
end

[M.Mean, ids, states, M.Unit] = sbtab_load_quantity_data(data_file, [], quantity_type, id_column, element_IDs, mean_columns,1);

if length(std_columns),
  M.Std = sbtab_load_quantity_data(data_file, [], quantity_type, id_column, element_IDs, std_columns,1);
else
  M.Std = [];
end
