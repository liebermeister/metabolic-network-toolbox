function ind_water = network_find_water(network)
  
% find metabolite index for metabolite 'water' based on possible annotations 
  
  ind_1 = label_names({'water'},network.metabolites);
  ind_2 = label_names({'H2O'},network.metabolites);
  if isfield(network,'Compound_Identifiers'),
    ind_3 = label_names('C00001',network.Compound_Identifiers);
    ind_4 = label_names('kegg:C00001',network.Compound_Identifiers);
  else
    ind_3=[]; ind_4=[];
  end
  if isfield(network,'metabolite_KEGGID'),
    ind_5 = find(strcmp('C00001',network.metabolite_KEGGID));
    ind_6 = find(strcmp('kegg:C00001',network.metabolite_KEGGID));
  else
    ind_5=0; ind_6=0;
  end
  ind = unique([ind_1, ind_2, ind_3, ind_4, ind_5, ind_6]);
  ind_water = ind(find(ind~=0));