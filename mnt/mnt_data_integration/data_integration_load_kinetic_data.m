function kinetic_data = data_integration_load_kinetic_data(data_quantities, quantity_info, network, file_kinetic_data, use_sbml_ids, use_kegg_ids, flag_invent_std, verbose)

% kinetic_data = data_integration_load_kinetic_data(data_quantities, quantity_info, network, file_kinetic_data, use_sbml_ids, use_kegg_ids, flag_invent_std, verbose)
%
% construct data structure with kinetic quantities relevant for a model
% for display: use parameter_balancing_kinetic_data_show(kinetic_data);
%
% data are read from file 'file_kinetic_data' and 
%
% the data structure 'kinetic_data' contains a field for every quantity type listed in 'data_quantities'
%
% fill in values from "kinetic_data_table"
% take logarithms where necessary
% values are NOT averaged yet (FIX!!!)

eval(default('quantity_info','[]','use_sbml_ids','1','use_kegg_ids','1','flag_invent_std','1','verbose','1'));

if isempty(quantity_info), quantity_info = data_integration_load_quantity_info; end

if iscell(file_kinetic_data),
  %% several files; read them one after the other, each time
  %% overriding information from previous files.
  kinetic_data = data_integration_load_kinetic_data(data_quantities, quantity_info, network, file_kinetic_data{1}, use_sbml_ids, use_kegg_ids);
  fn = fieldnames(kinetic_data);
  for it = 2:length(file_kinetic_data),
    my_kinetic_data = data_integration_load_kinetic_data(data_quantities, quantity_info, network, file_kinetic_data{it}, use_sbml_ids, use_kegg_ids);
    for it2 = 1:length(fn)
      ii = isfinite(my_kinetic_data.(fn{it2}).median);
      kinetic_data.(fn{it2}).median(ii) = my_kinetic_data.(fn{it2}).median(ii);
      ii = isfinite(my_kinetic_data.(fn{it2}).mean);
      kinetic_data.(fn{it2}).mean(ii) = my_kinetic_data.(fn{it2}).mean(ii);
      ii = isfinite(my_kinetic_data.(fn{it2}).std);
      kinetic_data.(fn{it2}).std(ii) = my_kinetic_data.(fn{it2}).std(ii);
      ii = isfinite(my_kinetic_data.(fn{it2}).lower);
      kinetic_data.(fn{it2}).lower(ii) = my_kinetic_data.(fn{it2}).lower(ii);
      ii = isfinite(my_kinetic_data.(fn{it2}).upper);
      kinetic_data.(fn{it2}).upper(ii) = my_kinetic_data.(fn{it2}).upper(ii);       

      if isfield(my_kinetic_data.(fn{it2}),'mean_ln'),
        ii = isfinite(my_kinetic_data.(fn{it2}).mean_ln);
        kinetic_data.(fn{it2}).mean_ln(ii) = my_kinetic_data.(fn{it2}).mean_ln(ii);
        ii = isfinite(my_kinetic_data.(fn{it2}).std_ln);
        kinetic_data.(fn{it2}).std_ln(ii) = my_kinetic_data.(fn{it2}).std_ln(ii);
        ii = isfinite(my_kinetic_data.(fn{it2}).lower_ln);
        kinetic_data.(fn{it2}).lower_ln(ii) = my_kinetic_data.(fn{it2}).lower_ln(ii);
        ii = isfinite(my_kinetic_data.(fn{it2}).upper);
        kinetic_data.(fn{it2}).upper_ln(ii) = my_kinetic_data.(fn{it2}).upper_ln(ii); 
      end
    end
  end
  return
end

if use_sbml_ids * isfield(network,'sbml_id_species'), 
  metabolites = network.sbml_id_species;
else,
  metabolites = network.metabolites;
end

if use_sbml_ids * isfield(network,'sbml_id_reaction'), 
  reactions   = network.sbml_id_reaction;
else,
  reactions   = network.actions;
end

% ------------------------------------------

if length(file_kinetic_data),
  kinetic_data_sbtab = sbtab_load_table(file_kinetic_data);
  kinetic_data_table = kinetic_data_sbtab.column.column;
else
  kinetic_data_table = struct('QuantityType',[],'SBMLSpeciesID',[]);
end

% ------------------------------------------

if use_kegg_ids, 
  if ~isfield(network,'metabolite_KEGGID'), 
    use_kegg_ids = 0;
    warning('Kegg IDs missing in network; matching elements by SBML ids');
  elseif ~isfield(kinetic_data_table,'Compound_MiriamID__urn_miriam_kegg_compound'), 
    use_kegg_ids= 0; 
    warning('Kegg IDs missing in data set; matching elements by SBML ids');
  else,
    metabolites = network.metabolite_KEGGID;
  end

  if isfield(network,'reaction_KEGGID'), 
    reactions = network.reaction_KEGGID;
  elseif  isfield(network,'MiriamID__urn_miriam_kegg_reaction'), 
    reactions = network.MiriamID__urn_miriam_kegg_reaction;
  else
    use_kegg_ids = 0; 
    warning('Kegg IDs missing in network; matching elements by SBML ids');
  end
end

% ------------------------------------------

[nm,nr] = size(network.N);

kinetic_data = struct;

for it = 1:length(data_quantities),
  
  flag_is_mean = 0;
  
  ind                     = find(strcmp(data_quantities{it},quantity_info.QuantityType));
  if isempty(ind), error(sprintf('Unknown quantity "%s"',data_quantities{it})); end
  symbol                  = quantity_info.Symbol{ind};
  scaling                 = quantity_info.Scaling{ind};
  related_element         = quantity_info.RelatedElement{ind};
  errstd                  = cell_string2num(quantity_info.ErrorStd(ind));
  quantity_entry          = struct;
  quantity_entry.scaling  = scaling;
  
  % construct empty vectors/matrices 

  switch quantity_info.RelatedElement{ind},
    case 'Species',          ss = [nm,1];
    case 'Reaction',         ss = [nr,1];
    case 'Reaction/Species', ss = [nr,nm];
    case 'None',             ss = [1];
  end
  
  quantity_entry.median = nan * zeros(ss);
  quantity_entry.mean   = nan * zeros(ss);
  quantity_entry.std    = nan * zeros(ss);
  quantity_entry.lower  = nan * zeros(ss);
  quantity_entry.upper  = nan * zeros(ss);

  % extract relevant kinetic data
  ind = find( strcmp(kinetic_data_table.QuantityType,data_quantities{it}));

  if isfield(kinetic_data_table,'Mean'),
    my_mean   = cell_string2num(kinetic_data_table.Mean(ind));
    my_median = nan * my_mean;
    flag_is_mean = 1;
  elseif isfield(kinetic_data_table,'Median'),
    my_median = cell_string2num(kinetic_data_table.Median(ind));
    my_mean   = nan * my_median;
  elseif isfield(kinetic_data_table,'Value'),
    my_median = cell_string2num(kinetic_data_table.Value(ind));
    my_mean   = nan * my_median;
  else 
    warning('No column with values found');
    my_median = [];
  end

  my_std    = nan* ones(length(ind),1);    
  my_std_ln = nan* ones(length(ind),1);    

   switch scaling,
    
    case 'Original',
      %% mean and median are the same
      if flag_is_mean,
        my_median = my_mean;
      else,
        my_mean = my_median;
      end
    
    case 'Logarithmic',
      my_mean_ln = nan * my_mean;
      my_std_ln  = nan * my_std;

  end

  if isfield(kinetic_data_table,'Std'),
    switch scaling,
      case 'Original',
        my_std    = cell_string2num(kinetic_data_table.Std(ind));
      case 'Logarithmic',
        my_std_ln  = log(10) * cell_string2num(kinetic_data_table.Std(ind));
    end
  end
    
  if flag_invent_std,
    if verbose, 
      display('Inventing standard deviations'); 
    end
    if length(my_mean),
    switch scaling,
      case 'Original',
        indices = find( [isfinite(my_mean) + isfinite(my_median)] .* [~isfinite(my_std)]);
        my_std(indices) = errstd;
      case 'Logarithmic',
        %% PRELIMINARY SOLUTION. FIX!
        indices = find( [isfinite(my_mean) + isfinite(my_median)] .* [~isfinite(my_std_ln)]);
        my_std_ln(indices)   = log(1+errstd);
    end
    end    
  end
  
  if isfield(kinetic_data_table,'LowerBound'),
    my_lower = cell_string2num(kinetic_data_table.LowerBound(ind));
  else,
    my_lower = nan * my_median;
  end

  if isfield(kinetic_data_table,'UpperBound'),
    my_upper = cell_string2num(kinetic_data_table.UpperBound(ind));
  else
    my_upper = nan * my_median;
  end

  switch related_element,
    case 'Species',   
      if use_kegg_ids, 
        rindices = label_names(kinetic_data_table.Compound_MiriamID__urn_miriam_kegg_compound(ind),metabolites,'multiple');
      else,
        rindices = label_names(kinetic_data_table.SBMLSpeciesID(ind),metabolites,'multiple');
      end
      cindices = repmat({1},length(rindices),1);

    case 'Reaction',  

      if use_kegg_ids, 
        if isfield(kinetic_data_table,'Reaction_MiriamID__urn_miriam_kegg_reaction'),
          rindices = label_names(kinetic_data_table.Reaction_MiriamID__urn_miriam_kegg_reaction(ind),reactions,'multiple');
          cindices = repmat({1},length(rindices),1);      
        else
          warning('IDs cannot be matched');
          rindices = repmat({},length(ind));
          cindices = repmat({},length(rindices),1);
        end
      elseif isfield(kinetic_data_table,'SBMLReactionID'),
        rindices = label_names(kinetic_data_table.SBMLReactionID(ind),reactions,'multiple');
        cindices = repmat({1},length(rindices),1);
      else
        warning('IDs cannot be matched');
        rindices = repmat({},length(ind));
        cindices = repmat({},length(rindices),1);
      end
  
    case 'Reaction/Species',   
      
      if use_kegg_ids, 
        if isfield(kinetic_data_table,'Reaction_MiriamID__urn_miriam_kegg_reaction'),
          rindices = label_names(kinetic_data_table.Reaction_MiriamID__urn_miriam_kegg_reaction(ind),reactions,'multiple');
          cindices = label_names(kinetic_data_table.Compound_MiriamID__urn_miriam_kegg_compound(ind),metabolites,'multiple');
        else
          warning('IDs cannot be matched');
          rindices = repmat({},length(ind));
          cindices = repmat({},length(rindices),1);
        end
      elseif isfield(kinetic_data_table,'SBMLReactionID'),
        rindices = label_names(kinetic_data_table.SBMLReactionID(ind),reactions,'multiple');
        cindices = label_names(kinetic_data_table.SBMLSpeciesID(ind),metabolites,'multiple');
      else
        warning('IDs cannot be matched');
        rindices = repmat({},length(ind));
        cindices = repmat({},length(rindices),1);
      end

    case 'None',   
      rindices = {};
      cindices = {};

  end
  
  ind_ok = find([cellfun('length',rindices)~=0].*[cellfun('length',cindices)~=0]);
  
  switch scaling,
    case 'Logarithmic',
      quantity_entry.mean_ln = nan * quantity_entry.mean;
      quantity_entry.std_ln  = nan * quantity_entry.mean;
  end
  
  %% this doesn't work for duplicate values: only the last value is accepted
  %% also doesn't work for values that have to appear in several places in the model
  
  for itt = 1:length(ind_ok),
    my_ind = ind_ok(itt);
    quantity_entry.mean(rindices{my_ind},cindices{my_ind})   = my_mean(my_ind);
    quantity_entry.median(rindices{my_ind},cindices{my_ind}) = my_median(my_ind);
    quantity_entry.std(rindices{my_ind}, cindices{my_ind})   = my_std(my_ind);
    quantity_entry.lower(rindices{my_ind},cindices{my_ind})  = my_lower(my_ind);
    quantity_entry.upper(rindices{my_ind}, cindices{my_ind}) = my_upper(my_ind);
    switch scaling,
      case 'Logarithmic',
        quantity_entry.mean_ln(rindices{my_ind}, cindices{my_ind}) = my_mean_ln(my_ind);
        quantity_entry.std_ln(rindices{my_ind}, cindices{my_ind})  = my_std_ln(my_ind);
    end
  end
  
  switch scaling,    
    case 'Logarithmic',
      if flag_is_mean,
        [quantity_entry.mean_ln, quantity_entry.std_ln] = lognormal_normal2log(quantity_entry.mean, quantity_entry.std);
        quantity_entry.median = exp(quantity_entry.mean_ln);
      else,
        quantity_entry.mean_ln = log(quantity_entry.median);
        ii = find(isfinite(quantity_entry.std_ln));
        [quantity_entry.mean(ii), quantity_entry.std(ii)] = lognormal_log2normal(quantity_entry.mean_ln(ii), quantity_entry.std_ln(ii));
      end
      quantity_entry.lower_ln = log(quantity_entry.lower);
      quantity_entry.upper_ln = log(quantity_entry.upper);
  end
  
  kinetic_data.(symbol) = quantity_entry; 
  
end
