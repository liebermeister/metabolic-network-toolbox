function kinetic_data = data_integration_load_kinetic_data(data_quantities, parameter_prior, network, data_file, use_sbml_ids, use_kegg_ids, flag_invent_std, verbose, filter_column, filter_entry, reaction_column_name, compound_column_name)

% kinetic_data = data_integration_load_kinetic_data(data_quantities, parameter_prior, network, data_file, use_sbml_ids, use_kegg_ids, flag_invent_std, verbose)
%
% Read model-related biochemical data from an SBtab file 
%
% Data are read from file 'data_file' 
%  - if no filename is provided, an empty data structure is generated)
%  - if 'data_file' contains a list of files, they are all read
%
% Output: 
%   Data structure 'kinetic_data' of kinetic quantities to be used in a model
%   It contains a field for every quantity type listed in 'data_quantities'
%   For displaying this data structure, use 'parameter_balancing_kinetic_data_show(kinetic_data);'

% fill in values from "kinetic_data_table"
% take logarithms where necessary
% values are NOT averaged but OVERWRITTEN!

if ~exist('sbtab_version','file'),
  error('For this function, the SBtab toolbox must be installed');
end

eval(default('data_quantities', '[]', ...
	     'parameter_prior','[]',...
	     'data_file', '[]', ...
             'use_sbml_ids','0',...
             'use_kegg_ids','0',...
             'flag_invent_std','1',...
             'verbose','0',...
             'filter_column','[]',...
             'filter_entry','[]', ...
             'reaction_column_name', '[]', ...
             'compound_column_name', '[]' ...
             ));

if isempty(data_quantities), 
data_quantities = {'standard Gibbs energy of reaction', ...
                   'standard chemical potential',...
                   'Michaelis constant',...
                   'activation constant', ...
                   'inhibitory constant',...
                   'equilibrium constant', ...
                   'substrate catalytic rate constant', ...
                   'product catalytic rate constant'};
end

if isempty(parameter_prior), 
  parameter_prior = biochemical_parameter_prior; 
end

if isempty(data_file),
  display('No data file provided. Creating empty data structure');
  verbose           = 0;
  data_file = [];
end


% ------------------------------------------
% if several data files are given, use all of them 

if iscell(data_file),

  %% several files; read them one after the other;
  %% each time information from previous files is overridden

  if isempty(reaction_column_name), reaction_column_name = repmat({[]},length(data_file),1); end
  if isempty(compound_column_name), compound_column_name = repmat({[]},length(data_file),1); end
  if prod(size(use_kegg_ids))==1, use_kegg_ids = use_kegg_ids * ones(length(data_file),1); end
  
  kinetic_data = data_integration_load_kinetic_data(data_quantities, parameter_prior, network, data_file{1}, use_sbml_ids, use_kegg_ids(1), flag_invent_std, verbose, filter_column, filter_entry,reaction_column_name{1},compound_column_name{1});
  fn = fieldnames(kinetic_data);

  for it = 2:length(data_file),
    my_kinetic_data = data_integration_load_kinetic_data(data_quantities, parameter_prior, network, data_file{it}, use_sbml_ids, use_kegg_ids(2), flag_invent_std, verbose, filter_column, filter_entry,reaction_column_name{it},compound_column_name{it});
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


% ------------------------------------------
% if filter criterion is given, put the rows that meet this criterion
% at the end, such that they will be used!

if length(data_file),
  kinetic_data_sbtab = sbtab_table_load(data_file);
  kinetic_data_table = sbtab_table_get_all_columns(kinetic_data_sbtab);
  if length(filter_column),
    if isfield(kinetic_data_table,filter_column),
      display(sprintf('File %s:\n  Preferring entries marked as %s in column %s', data_file, filter_entry, filter_column));
      ind_good_entries = find(strcmp(kinetic_data_table.(filter_column),filter_entry));
      ind_bad_entries  = find(strcmp(kinetic_data_table.(filter_column),filter_entry)==0);
      my_fn = fieldnames(kinetic_data_table);
      for it = 1:length(my_fn),
        kinetic_data_table.(my_fn{it}) = kinetic_data_table.(my_fn{it})([ind_bad_entries; ind_good_entries]);
      end
    end
  else
    if length(filter_column),
      display(sprintf('No column %s found in file %s - no preference for entries', filter_column, data_file));
    end
  end
else
  % create empty table
  kinetic_data_table = struct('QuantityType',[],'Compound_SBML_species_id',[]);
end


% ------------------------------------------

if use_sbml_ids,
  if ~isfield(kinetic_data_table,'Reaction_SBML_reaction_id'),
    if ~isfield(kinetic_data_table,'Compound_SBML_species_id'),
      warning('No SBML IDs found in data file');
      use_sbml_ids = 0; 
    end
  end
end

if use_kegg_ids,
  if ~isfield(kinetic_data_table,'Reaction_Identifiers_kegg_reaction'),
    if ~isfield(kinetic_data_table,'Compound_Identifiers_kegg_compound'),
      warning('No KEGG IDs found in data file');
      use_kegg_ids = 0; 
    end
  end
end

if isempty(reaction_column_name),
  if use_sbml_ids,
    reaction_column_name = 'Reaction_SBML_reaction_id';
  elseif use_kegg_ids,
    reaction_column_name = 'Reaction_Identifiers_kegg_reaction';
  else
    reaction_column_name = 'Reaction';
  end
end

if isempty(compound_column_name),
  if use_sbml_ids,
    reaction_column_name = 'Compound_SBML_species_id';
  elseif use_kegg_ids,
    compound_column_name =  'Compound_Identifiers_kegg_compound';
  else
    compound_column_name =  'Compound';
  end
end


% ------------------------------------------

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

if use_kegg_ids,
  
  if ~isfield(network,'metabolite_KEGGID'), 
    if isfield(network,'Identifiers_kegg_compound'), 
      network.metabolite_KEGGID = network.Identifiers_kegg_compound;
    end
  end
      
  if ~isfield(network,'reaction_KEGGID'), 
    if isfield(network,'Identifiers_kegg_reaction'), 
      network.reaction_KEGGID = network.Identifiers_kegg_reaction;
    end
  end
      
  if ~isfield(network,'metabolite_KEGGID'), 
    % use_kegg_ids = 0;
    if verbose, 
      warning('Metabolite Kegg IDs missing in network');%; matching elements by SBML ids');
    end
  elseif ~isfield(kinetic_data_table,compound_column_name), 
    %use_kegg_ids= 0;
    if verbose, 
      warning(sprintf('%s: Compound Kegg IDs missing in data set',data_file));%; matching elements by SBML ids');
    end
  else,
    metabolites = network.metabolite_KEGGID;
  end

  if isfield(network,'reaction_KEGGID'), 
    reactions = network.reaction_KEGGID;
  elseif isfield(network,'Identifiers_kegg_reaction'), 
    reactions = network.Identifiers_kegg_reaction;
  else
    %use_kegg_ids = 0;
    if verbose, 
      warning('Reaction Kegg IDs missing in network');%; matching elements by SBML ids');
    end
  end
end

% ------------------------------------------

[nm,nr] = size(network.N);

kinetic_data = struct;

for it = 1:length(data_quantities),
  
  flag_is_mean = 0;
  
  ind                     = find(strcmp(data_quantities{it},parameter_prior.QuantityType));
  if isempty(ind), error(sprintf('Unknown quantity "%s"',data_quantities{it})); end
  symbol                  = parameter_prior.Symbol{ind};
  scaling                 = parameter_prior.Scaling{ind};
  related_element         = parameter_prior.Element{ind};
  errstd                  = cell_string2num(parameter_prior.DataStd(ind));
  quantity_entry          = struct;
  quantity_entry.scaling  = scaling;
  
  % construct empty vectors/matrices 

  switch parameter_prior.Element{ind},
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

  if isfield(kinetic_data_table,'Median'),
    my_median = cell_string2num(kinetic_data_table.Median(ind));
    my_mean   = nan * my_median;
  elseif isfield(kinetic_data_table,'Mean'),
    my_mean   = cell_string2num(kinetic_data_table.Mean(ind));
    my_median = nan * my_mean;
    flag_is_mean = 1;
  elseif isfield(kinetic_data_table,'Value'),
    my_median = cell_string2num(kinetic_data_table.Value(ind));
    my_mean   = nan * my_median;
  else 
    if verbose, warning('No column with values found'); end
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
      display(sprintf('  Quantity %s: Inventing standard deviations',symbol)); 
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
    error('Deprecated column LowerBound found; please fix the data file and use Min and Max instead.');
  end

  if isfield(kinetic_data_table,'Min'),
    my_lower = cell_string2num(kinetic_data_table.Min(ind));
  else,
    my_lower = nan * my_median;
  end

  if isfield(kinetic_data_table,'Max'),
    my_upper = cell_string2num(kinetic_data_table.Max(ind));
  else
    my_upper = nan * my_median;
  end

  switch related_element,
    
    case 'Species',  
      if use_kegg_ids,
        if isfield(kinetic_data_table,compound_column_name),
          rindices = label_names(kinetic_data_table.(compound_column_name)(ind),metabolites,'multiple');
        else
          rindices = {};
        end
      else,
        rindices = label_names(kinetic_data_table.Compound_SBML_species_id(ind),metabolites,'multiple');
      end
      cindices = repmat({1},length(rindices),1);
      
    case 'Reaction',

      if use_kegg_ids,
        
        if isfield(kinetic_data_table,reaction_column_name),
          rindices = label_names(kinetic_data_table.(reaction_column_name)(ind),reactions,'multiple');
        else
          if verbose,
            if isfield(kinetic_data, symbol),
              warning(sprintf('File %s: Looking for %s data: Reaction IDs cannot be matched',data_file,symbol));
            end
          end
          rindices = repmat({1},length(ind));
        end

        if isfield(kinetic_data_table,compound_column_name),          
          cindices = repmat({1},length(rindices),1);      
        else
          if verbose, 
            if isfield(kinetic_data, symbol),
              warning(sprintf('File %s: Looking for %s data: Compound IDs cannot be matched',data_file,symbol));
            end
          end
          cindices = repmat({1},length(rindices),1);
        end
        
      elseif isfield(kinetic_data_table,'Reaction_SBML_reaction_id'),
        rindices = label_names(kinetic_data_table.Reaction_SBML_reaction_id(ind),reactions,'multiple');
        cindices = repmat({1},length(rindices),1);
      
      else
        if verbose, 
          warning(sprintf('%s: %s: Reaction IDs cannot be matched',data_file,symbol));
        end
        rindices = repmat({1},length(ind));
        cindices = repmat({1},length(rindices),1);
      end
      
    case 'Reaction/Species',   
      if use_kegg_ids,
        if isfield(kinetic_data_table,reaction_column_name) * ...
              isfield(kinetic_data_table,compound_column_name),
          rindices = label_names(kinetic_data_table.(reaction_column_name)(ind),reactions,'multiple');
          cindices = label_names(kinetic_data_table.(compound_column_name)(ind),metabolites,'multiple');
        else
          if verbose, 
            warning(sprintf('%s: %s: Joint Compound and Reaction IDs cannot be matched',data_file,symbol));
          end         
          rindices = repmat({1},length(ind),1);
          cindices = repmat({1},length(rindices),1);
        end
      elseif isfield(kinetic_data_table,'Reaction_SBML_reaction_id'),
        rindices = label_names(kinetic_data_table.Reaction_SBML_reaction_id(ind),reactions,'multiple');
        cindices = label_names(kinetic_data_table.Compound_SBML_species_id(ind),metabolites,'multiple');
      else
        if verbose,
          warning(sprintf('%s: %s: Joint Compound and Reaction IDs cannot be matched',data_file,symbol));
        end
        rindices = repmat({1},length(ind));
        cindices = repmat({1},length(rindices),1);
      end

    case 'None',   
      rindices = {};
      cindices = {};

  end 

  if isempty(rindices), rindices ={}; end
  if isempty(cindices), cindices ={}; end

  %% -----------------------------------------------------------

  ind_ok = find([cellfun('length',rindices)~=0].*[cellfun('length',cindices)~=0]);

  % symbol
  % rindices
  % kinetic_data_table.Value(ind(ind_ok))
  % kinetic_data_table.Reaction_Identifiers_kegg_reaction(ind(ind_ok))
  
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


% --------------------------------------------------------------------------------------
% make sure important fields in kinetic_data are considered


if sum(isfinite(kinetic_data.Keq.lower_ln )) == 0, 
  kinetic_data.Keq.lower_ln = log(kinetic_data.Keq.lower);
end

if sum(isfinite(kinetic_data.Keq.upper_ln )) == 0, 
  kinetic_data.Keq.upper_ln = log(kinetic_data.Keq.upper);
end

if sum(isfinite(kinetic_data.Kcatf.lower_ln )) == 0, 
  kinetic_data.Kcatf.lower_ln = log(kinetic_data.Kcatf.lower);
end

if sum(isfinite(kinetic_data.Kcatf.upper_ln )) == 0, 
  kinetic_data.Kcatf.upper_ln = log(kinetic_data.Kcatf.upper);
end

if sum(isfinite(kinetic_data.Kcatr.lower_ln )) == 0, 
  kinetic_data.Kcatr.lower_ln = log(kinetic_data.Kcatr.lower);
end

if sum(isfinite( kinetic_data.Kcatr.upper_ln )) == 0, 
  kinetic_data.Kcatr.upper_ln = log(kinetic_data.Kcatr.upper);
end

