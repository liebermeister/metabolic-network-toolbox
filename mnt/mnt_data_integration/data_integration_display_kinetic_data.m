function T = data_integration_display_kinetic_data(kinetic_data,network,flag_display)

% T = data_integration_display_kinetic_data(kinetic_data,network,flag_display)

eval(default('flag_display','1'));

fn = fieldnames(kinetic_data);

for it = 1:length(fn),

  xx = kinetic_data.(fn{it});
  
  %ii = find(isfinite(xx.median) + isfinite(xx.mean) + isfinite(xx.std));
  ii = find(isfinite(xx.median) + isfinite(xx.mean) + isfinite(xx.std) + isfinite(xx.lower) + isfinite(xx.upper));
  nii = length(ii);
  
  T{it} = [{'!Median','!Mean','!Std','!Lower','!Upper'}; ...
           num2cell([xx.median(ii), ...
                     xx.mean(ii), ...
                     xx.std(ii), ...
                     xx.lower(ii), ...
                     xx.upper(ii)])];
  
    quantity_info   = data_integration_load_quantity_info;
    ind = find(strcmp(fn{it},quantity_info.Symbol));
    related_element      = quantity_info.RelatedElement{ind};
    quantity_type        = quantity_info.QuantityType{ind};
    quantity_type_column = [{'!QuantityType'}; repmat({quantity_type},size(T{it},1)-1,1)];
    unit                 = quantity_info.Unit{ind};
    unit_column          = [{'!Unit'}; repmat({unit},size(T{it},1)-1,1)];
    
    [nm,nr] = size(network.N);
    
    species_kegg_column  = [];
    reaction_kegg_column = [];
    
    switch related_element,
      
      case 'Species',
        species_column       = [ {'!Compound'};  network.metabolites(ii)];
        if isfield(network, 'metabolite_KEGGID'),
          species_kegg_column = [{'!CompoundID:urn.miriam.kegg.compound'}; network.metabolite_KEGGID(ii)];
        end
        reaction_column      = [{'!Reaction'}; repmat({''},nii,1)];
        if isfield(network, 'MiriamID_urn_miriam_kegg_reaction'),
          reaction_kegg_column = [{'!ReactionID:urn.miriam.kegg.reaction'}; repmat({''},nii,1)];
        end
        
      case 'Reaction', 
        species_column  = [{'!Compound'}; repmat({''},nii,1)];
        reaction_column = [{'!Reaction'}; network.actions(ii)];
        if isfield(network, 'metabolite_KEGGID'),
          species_kegg_column = [{'!CompoundID:urn.miriam.kegg.compound'}; repmat({''},nii,1)];
        end
        if isfield(network, 'MiriamID_urn_miriam_kegg_reaction'),
          reaction_kegg_column = [{'!ReactionID:urn.miriam.kegg.reaction'}; network.MiriamID_urn_miriam_kegg_reaction];
        end
    
      case 'Reaction/Species',
        [ij,ik] = ind2sub(size(network.N),ii);
        species_column       = [ {'!Compound'};  network.metabolites(ij)];
        if isfield(network, 'metabolite_KEGGID'),
          species_kegg_column = [{'!CompoundID:urn.miriam.kegg.compound'}; network.metabolite_KEGGID(ij)];
        end
        reaction_column = [{'!Reaction'}; network.actions(ik)];
        if isfield(network, 'MiriamID_urn_miriam_kegg_reaction'),
          reaction_kegg_column = [{'!ReactionID:urn.miriam.kegg.reaction'}; network.MiriamID_urn_miriam_kegg_reaction(ik)];
        end
        
      case 'None',
        species_kegg_column  = [];
        reaction_kegg_column = [];

    end
    if nii,
      T{it} = [ quantity_type_column, species_column, reaction_column, unit_column, T{it}, species_kegg_column, reaction_kegg_column];
    end
  if flag_display, 
    table(T{it},0)
  end
end
