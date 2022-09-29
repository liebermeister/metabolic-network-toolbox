function kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, data_file, options)

% kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, data_file, options)
%
% Read biochemical (kinetic or thermodynamic) data from SBtab file and map them to network model
%
% Data are read from file 'data_file' 
%  - if no filename is provided, an empty data structure is generated)
%  - if 'data_file' contains a list of files, they are all used
%
% Output: 
%   Data structure 'kinetic_data' of kinetic quantities (to be used in a model)
%   It contains a field for every quantity type listed in 'data_quantities'
%
% Options
%   options.use_sbml_ids
%   options.use_kegg_ids
%   options.flag_invent_std
%   options.verbose
%   options.filter_column
%   options.filter_entry
%   options.reaction_column_name
%   options.compound_column_name
%   options.use_python_version_defaults
%   options.expected_unit
%   options.enforce_ranges
%   options.table_id (optional, table ID if the input file contains several tables)
%
% Other functions for kinetic data structure:
%   Display data values:             parameter_balancing_kinetic_data_show.m
%   Save data to file  :             kinetic_data_save.m
%   Insert bounds and pseudo values: kinetic_data_complete.m
%   Select data for a submodel:      kinetic_data_subselect.m
 
% values are NOT averaged but OVERWRITTEN!
% fill in values from "kinetic_data_table"
% take logarithms where necessary
  
if ~exist('sbtab_version','file'), error('For this function, the SBtab toolbox must be installed'); end

eval(default('data_quantities', '[]', 'parameter_prior','[]', 'data_file', '[]', 'options','struct'));

options_default = struct('use_sbml_ids',0, ...
                         'use_kegg_ids',0, ...
                         'flag_invent_std',1, ...
                         'verbose',0, ...
                         'filter_column',[], ...
                         'filter_entry',[], ...
                         'reaction_column_name', [], ...
                         'compound_column_name', [], ...
                         'use_python_version_defaults', 0, ...
                         'expected_unit',[],...
                         'enforce_ranges',1,...
                         'table_id',[]);

options = join_struct(options_default,options);

verbose = options.verbose;

if isempty(parameter_prior), 
  parameter_prior = parameter_balancing_prior; 
end

if isempty(data_quantities), 
  data_quantities = {'standard Gibbs free energy of reaction', ...
                     'standard chemical potential',...
                     'Michaelis constant',...
                     'activation constant', ...
                     'inhibitory constant',...
                     'equilibrium constant', ...
                     'substrate catalytic rate constant', ...
                     'product catalytic rate constant'};
end

flag_construct_empty_data_structure = 0;

if isempty(data_file),
  if verbose,
    display('o No data file provided. Creating empty data structure');
  end
  flag_construct_empty_data_structure = 1;
  verbose = 0;
  data_file       = [];
elseif isstr(data_file),
  if verbose, 
    display(sprintf('o Collecting data from file %s', data_file));
  end
end


% ------------------------------------------
% if several data files are given, use all of them 

if iscell(data_file),

  %% several files; read them one after the other;
  %% each time information from previous files is overridden

  if isempty(options.reaction_column_name), options.reaction_column_name = repmat({[]},length(data_file),1); end
  if isempty(options.compound_column_name), options.compound_column_name = repmat({[]},length(data_file),1); end
  if prod(size(options.use_kegg_ids))==1,   options.use_kegg_ids = options.use_kegg_ids * ones(length(data_file),1); end
  
  my_opt = options;
  if isstr(options.reaction_column_name),
    my_opt.reaction_column_name = options.reaction_column_name;
  else
    my_opt.reaction_column_name = options.reaction_column_name{1};
  end
  my_opt.compound_column_name = options.compound_column_name{1};
  my_opt.table_id = options.table_id;
  kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, data_file{1}, my_opt);
  fn = fieldnames(kinetic_data);

  for it = 2:length(data_file),
    my_opt = options;
    my_opt.use_kegg_ids = options.use_kegg_ids(2);
    if isstr(options.reaction_column_name),
      my_opt.reaction_column_name = options.reaction_column_name;
    else
      my_opt.reaction_column_name = options.reaction_column_name{1};
    end
    my_opt.compound_column_name = options.compound_column_name{it};
    my_kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, data_file{it}, my_opt);
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
  if isstr(data_file),
    kinetic_data_sbtab = sbtab_table_load(data_file,[],[],options.table_id);
  else
    kinetic_data_sbtab = data_file; data_file='[...]';
  end
  kinetic_data_table = sbtab_table_get_all_columns(kinetic_data_sbtab);
  
  ll = find(strcmp('forward catalytic rate constant',kinetic_data_table.QuantityType));
  kinetic_data_table.QuantityType(ll) = repmat({'substrate catalytic rate constant'},length(ll),1);
  ll = find(strcmp('reverse catalytic rate constant',kinetic_data_table.QuantityType));
  kinetic_data_table.QuantityType(ll) = repmat({'product catalytic rate constant'},length(ll),1);
                      
  if ~isfield(kinetic_data_table,'Unit'),
    error(sprintf('Column "Unit" missing in data file %s',data_file));
  end
  if length(options.filter_column),
    if isfield(kinetic_data_table,options.filter_column),
        if verbose,
      display(sprintf('File %s:\n  Preferring entries marked as %s in column %s', data_file, options.filter_entry, options.filter_column));
      end
      ind_good_entries = find(strcmp(kinetic_data_table.(options.filter_column),options.filter_entry));
      ind_bad_entries  = find(strcmp(kinetic_data_table.(options.filter_column),options.filter_entry)==0);
      my_fn = fieldnames(kinetic_data_table);
      for it = 1:length(my_fn),
        kinetic_data_table.(my_fn{it}) = kinetic_data_table.(my_fn{it})([ind_bad_entries; ind_good_entries]);
      end
    end
  else
    if length(options.filter_column),
      if verbose,
        display(sprintf('No column %s found in file %s - no preference for entries', options.filter_column, data_file));
      end
    end
  end
else
  % create empty table
  kinetic_data_table = struct('QuantityType',[],'Compound_SBML_species_id',[]);
  kinetic_data_sbtab = [];
end


% ------------------------------------------
% choose reaction and compound ID columns to use

kegg_ids_found  = isfield(kinetic_data_table,'Reaction_Identifiers_kegg_reaction') + isfield(kinetic_data_table,'Compound_Identifiers_kegg_compound');
sbml_ids_found  = isfield(kinetic_data_table,'Reaction_SBML_reaction_id') + isfield(kinetic_data_table,'Compound_SBML_species_id');
sbtab_ids_found = isfield(kinetic_data_table,'Reaction') + isfield(kinetic_data_table,'Compound');

use_IDs = 'none';
if sbtab_ids_found;
  use_IDs = 'sbtab';
end
if options.use_sbml_ids,
  if sbml_ids_found;
    use_IDs = 'sbml';
  else
    if verbose,
      display(sprintf('- (kinetic_data_load.m): No SBML IDs found in data file', data_file));
    end
  end
end
if options.use_kegg_ids,
  if kegg_ids_found, 
    use_IDs = 'kegg';
  else
    if verbose,
      display(sprintf('- (kinetic_data_load.m): No KEGG IDs found in data file', data_file));  
    end
  end
end

if strcmp(use_IDs,'none'),
  if sbml_ids_found, 
    use_IDs = 'sbml';
  end
  if sbtab_ids_found, 
    use_IDs = 'sbtab';
  end
end

if strcmp(use_IDs,'none'),
  error('I don not know what ID columns to use');
end

switch use_IDs,
  case 'kegg'
    if ~kegg_ids_found,
      options.use_kegg_ids = 0; 
    end
  case 'sbml',
    if ~sbml_ids_found,
      if verbose,
        display(sprintf('- (kinetic_data_load.m): No SBML IDs found in data file %s; trying to use SBtab IDs instead', data_file));
      end
      options.use_sbml_ids = 0;
      use_sbtab_ids = 1;
    end
  case 'sbtab',
    if ~sbtab_ids_found,
      if verbose,
        display(sprintf('- (kinetic_data_load.m): No SBML IDs found in data file %s; trying to use SBML IDs instead', data_file));
      end
      options.use_sbml_ids = 1;
    end
end

%if isempty(options.reaction_column_name),
  switch use_IDs, 
    case 'kegg',
      options.reaction_column_name = 'Reaction_Identifiers_kegg_reaction';
    case 'sbml',
      options.reaction_column_name = 'Reaction_SBML_reaction_id';
    case 'sbtab'
      options.reaction_column_name = 'Reaction';
  end
%end

%if isempty(options.compound_column_name),
  switch use_IDs, 
    case 'kegg',
      options.compound_column_name = 'Compound_Identifiers_kegg_compound';
    case 'sbml',
      options.compound_column_name = 'Compound_SBML_species_id';
    case 'sbtab'
      options.compound_column_name = 'Compound';
  end
%end

% ------------------------------------------

if options.use_sbml_ids * isfield(network,'sbml_id_species'), 
  metabolites = network.sbml_id_species;
else,
  metabolites = network.metabolites;
end

if options.use_sbml_ids * isfield(network,'sbml_id_reaction'), 
  reactions   = network.sbml_id_reaction;
else,
  reactions   = network.actions;
end


% ------------------------------------------

if options.use_kegg_ids,
  
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
    if verbose, display('- (kinetic_data_load.m): Metabolite KEGG IDs missing in network'); end
  elseif ~isfield(kinetic_data_table,options.compound_column_name), 
    if verbose, display(sprintf('- (kinetic_data_load.m): No compound KEGG IDs found')); end
  else,
    metabolites = network.metabolite_KEGGID;
  end

  if ~isfield(network,'reaction_KEGGID'), 
    if verbose, display('- (kinetic_data_load.m): Reaction KEGG IDs missing in network'); end
  elseif ~isfield(kinetic_data_table,options.reaction_column_name), 
    if verbose, display(sprintf('- (kinetic_data_load.m): No reaction KEGG IDs found')); end
  else
    reactions = network.reaction_KEGGID;
  end

end

% ------------------------------------------

[nm,nr] = size(network.N);

kinetic_data = struct;

for it = 1:length(data_quantities),
  flag_mean_is_given      = 0;
  ind                     = find(strcmp(data_quantities{it},parameter_prior.QuantityType));
  if isempty(ind), 
    error(sprintf('Unknown quantity "%s"',data_quantities{it})); 
  end
  symbol                  = parameter_prior.Symbol{ind};
  quantitytype            = parameter_prior.QuantityType{ind};
  scaling                 = parameter_prior.MathematicalType{ind};
  related_element         = parameter_prior.BiologicalElement{ind};
  default_unit            = parameter_prior.Unit{ind};
  if length(options.expected_unit),
    default_unit = options.expected_unit;
  end
  errstd      = cell_string2num(parameter_prior.DataStd(ind));
  errgeomstd  = cell_string2num(parameter_prior.DataGeometricStd(ind));
  allowed_min = eval(parameter_prior.LowerBound{ind});
  allowed_max = eval(parameter_prior.UpperBound{ind});

% construct empty vectors/matrices 

  require_compound_id = 0;
  require_reaction_id = 0;

  switch parameter_prior.BiologicalElement{ind},
    case 'Species',          ss = [nm,1];   require_compound_id=1;
    case 'Reaction',         ss = [nr,1];   require_reaction_id=1;
    case 'Reaction/Species', ss = [nr,nm];  require_compound_id=1; require_reaction_id=1;
    case 'None',             ss = [1];
  end
  
  if verbose,
    if ~isempty(data_file),
      if require_compound_id,
        if ~isfield(kinetic_data_table,options.compound_column_name),
          warning(sprintf('No compound ID column %s found in data file %s',options.compound_column_name, data_file));
        end
      end
      if require_reaction_id,
        if ~isfield(kinetic_data_table,options.reaction_column_name),
          warning(sprintf('No reaction ID column %s found in data file %s',options.reaction_column_name, data_file));
        end
      end
    end
  end
  
  %% extract relevant kinetic data
  
  %kinetic_data_table.QuantityType
  ind = find(strcmp(kinetic_data_table.QuantityType,data_quantities{it}));
    
  if ~isempty(data_file),
    %% omit data with non-standard units
    my_units = kinetic_data_table.Unit(ind);
    is_ok = double(strcmp(my_units,default_unit));
    if sum(is_ok==0), 
      if verbose,
      display(sprintf('- (kinetic_data_load.m): The data file contains "%s" values with unknown units. I ignore these values.',quantitytype)); 
      end
      data_file
      default_unit
      units_found = my_units
    end
    ind = ind(find(is_ok));
  end
  
  my_median  = nan * ones(length(ind),1);
  my_mean    = nan * ones(length(ind),1);
  my_std     = nan * ones(length(ind),1);
  my_mean_ln = nan * ones(length(ind),1);    
  my_std_ln  = nan * ones(length(ind),1);    
  my_lower   = nan * ones(length(ind),1);
  my_upper   = nan * ones(length(ind),1);
  
  flag_mean_is_given = 0;
  
  if isfield(kinetic_data_table,'Median'),
    my_median = cell_string2num(kinetic_data_table.Median(ind));
  elseif isfield(kinetic_data_table,'Value'),
    my_median = cell_string2num(kinetic_data_table.Value(ind));
  elseif isfield(kinetic_data_table,'Mode'),
    my_median = cell_string2num(kinetic_data_table.Mode(ind));
  elseif ~isfield(kinetic_data_table,'Mean'),
    if ~isfield(kinetic_data_table,'Min'),
      if ~isfield(kinetic_data_table,'Max'),
        if length(data_file),
          warning('- (kinetic_data_load.m): No column with data values found'); 
        end
      end
    end
  end
  
  if isfield(kinetic_data_table,'Mean'),
    flag_mean_is_given = 1;
    my_mean = cell_string2num(kinetic_data_table.Mean(ind));
    if prod(isnan(my_mean)),
      flag_mean_is_given = 0;
    end
  end

  if options.enforce_ranges,
    
    ind_mean_too_low    = find(my_mean   < allowed_min);
    ind_mean_too_high   = find(my_mean   > allowed_max);
    ind_median_too_low  = find(my_median < allowed_min);
    ind_median_too_high = find(my_median > allowed_max);
    
    if length(ind_mean_too_low) + length(ind_mean_too_high) + length(ind_median_too_low) + length(ind_median_too_high),
      if verbose,
        display(sprintf('- (kinetic_data_load.m): %s data values outside allowed range in file %s. \n  I replace these values by the minimal or maximal allowed values %f and %f', symbol, data_file,allowed_min, allowed_max));
        display([my_mean(ind_mean_too_low); my_mean(ind_mean_too_high); my_median(ind_median_too_low); my_median(ind_median_too_high)])
      end
      my_mean(ind_mean_too_low)      = allowed_min;
      my_mean(ind_mean_too_high)     = allowed_max;
      my_median(ind_median_too_low)  = allowed_min;
      my_median(ind_median_too_high) = allowed_max;
    end
  
  end
  
  switch scaling,
    case 'Multiplicative',
      ind_zero_values = find([my_mean==0]+[my_median==0]);
      if length(ind_zero_values),
        if verbose,
          display(sprintf('- (kinetic_data_load.m): The data for quantity "%s" contain zero values. I will ignore them', symbol));
        end
        my_mean(ind_zero_values)  = nan;
        my_median(ind_zero_values)= nan;
      end
  end

  if isfield(kinetic_data_table,'Std'),
    my_std = cell_string2num(kinetic_data_table.Std(ind));
  end
  ind_zero_values = find(my_std==0);
  if length(ind_zero_values),
    if verbose,
      display(sprintf('- (kinetic_data_load.m): The data for quantity "%s" contain zero standard deviations. I will ignore them', symbol));
    end
    my_std(ind_zero_values)=nan;
  end
  
  if options.flag_invent_std,
    indices_std_missing = find( [isfinite(my_mean) + isfinite(my_median)] .* [~isfinite(my_std)]);
    if verbose, 
      if length(indices_std_missing),
        if verbose,
          display(sprintf('    Quantity %s: Inventing standard deviations',symbol)); 
        end
      end
    end

    %if options.use_python_version_defaults,
      %% Possibility 1: complete missing data error by using default geometric standard deviation
      switch scaling,
        case 'Additive',
          my_std_guess = errstd * ones(size(my_median));
        case 'Multiplicative',
          if flag_mean_is_given,
            my_std_guess = my_mean * sqrt(exp(log(errgeomstd).^2)-1); % formula assuming log-normal distribution
          else
            [~,my_std_guess] = lognormal_log2normal(log(my_median),log(errgeomstd) * ones(size(my_median)));
          end
      end
    %else
    %  %% Possibility 2: complete missing data error by default standard deviation
    %  my_std_guess = errstd * ones(size(my_std));
    %end
    my_std(indices_std_missing) = my_std_guess(indices_std_missing);
  end

  switch scaling,
    case 'Additive',
      if flag_mean_is_given, my_median = my_mean;
      else,                  my_mean = my_median;
      end
      my_geomstd = nan * my_mean;
    case 'Multiplicative',
      if flag_mean_is_given,
        [my_mean_ln,my_std_ln] = lognormal_normal2log(my_mean,my_std);
        my_median(isfinite(my_mean_ln)) = exp(my_mean_ln(isfinite(my_mean_ln)));
        my_geomstd(isfinite(my_std_ln)) = exp(my_std_ln(isfinite(my_std_ln)));        
      else,
        %% THIS IS A FIX! Use predefined geometric standard deviation
        my_mean_ln = log(my_median);
        my_std_ln  = log(errgeomstd) * ones(size(my_mean_ln));
        [my_mean,my_std] = lognormal_log2normal(my_mean_ln,my_std_ln);
        my_geomstd(isfinite(my_std_ln)) = exp(my_std_ln(isfinite(my_std_ln)));        
      end
  end

  if isfield(kinetic_data_table,'Min'),
    my_lower = cell_string2num(kinetic_data_table.Min(ind));
  end
  
  if isfield(kinetic_data_table,'Max'),
    my_upper = cell_string2num(kinetic_data_table.Max(ind));
  end
  
  %% Determine reaction and species indices in the model
  
  rindices = [];
  cindices = [];

  if length(kinetic_data_table.QuantityType),

  switch related_element,

    case 'Species',
      if isfield(kinetic_data_table,options.compound_column_name),
        rindices = label_names(kinetic_data_table.(options.compound_column_name)(ind),metabolites);%,'multiple');
      else
        if verbose,
          display(sprintf('- (kinetic_data_load.m): Looking for %s data: Compound IDs cannot be matched',symbol));
        end
        rindices = repmat(0,length(ind),1);
        if verbose, 
          if isfield(kinetic_data, symbol),
            if verbose,
              display(sprintf('- (kinetic_data_load.m): Looking for %s data: Compound IDs cannot be matched',symbol));
            end
          end
        end
      end
      cindices = repmat(1,length(rindices),1);
      if find(rindices==0),
        if verbose,
          display(sprintf('- (kinetic_data_load.m): Unknown species ID for %s encountered in data file',symbol));
        end
        if verbose,
          unknown_species_ids = mytable(kinetic_data_table.(options.compound_column_name)(ind(find(rindices==0))),0)
        else
          if verbose,
            display('Set option "verbose" to see the species ID');
          end
        end
      end
      
    case 'Reaction',

      if isfield(kinetic_data_table,options.reaction_column_name),
        rindices = label_names(kinetic_data_table.(options.reaction_column_name)(ind),reactions);%,'multiple');
      else
        if verbose,
          display(sprintf('- (kinetic_data_load.m): Looking for %s data: Reaction IDs cannot be matched',symbol));
        end
        rindices = repmat(1,length(ind),1);
        if verbose,
          if isfield(kinetic_data, symbol),
            if verbose,
              display(sprintf('- (kinetic_data_load.m): Looking for %s data: Reaction IDs cannot be matched',symbol));
            end
          end
        end
      end
      cindices = repmat(1,length(rindices),1);      
      if find(rindices==0),
        if verbose,
          display(sprintf('- (kinetic_data_load.m): Unknown reaction ID for %s encountered in data file',symbol));
        end
        if verbose,
          unknown_reaction_ids = mytable(kinetic_data_table.(options.reaction_column_name)(ind(find(rindices==0))),0)
        end
      end
      
    case 'Reaction/Species',   
      if isfield(kinetic_data_table,options.reaction_column_name) * isfield(kinetic_data_table,options.compound_column_name),
        rindices = label_names(kinetic_data_table.(options.reaction_column_name)(ind),reactions);%,  'multiple');
        cindices = label_names(kinetic_data_table.(options.compound_column_name)(ind),metabolites);%,'multiple');
      else
        if verbose, 
          display(sprintf('    Quantity %s: Joint Compound and Reaction IDs cannot be matched',symbol)); 
        end         
        rindices = repmat(1,length(ind),1);
        cindices = repmat(1,length(rindices),1);
      end
      if find(rindices==0),
        if verbose,
          display(sprintf('- (kinetic_data_load.m): Unknown reaction ID for %s encountered in data file',symbol));
        end
        if verbose,
          unknown_reaction_ids = mytable(kinetic_data_table.(options.reaction_column_name)(ind(find(rindices==0))),0)
          species_ids = mytable(kinetic_data_table.(options.compound_column_name)(ind(find(rindices==0))),0)
        else
          if verbose,
            display('Set option "verbose" to see the reaction ID');
          end
        end
      end
      if find(cindices==0),
        if verbose,
          display(sprintf('- (kinetic_data_load.m): Unknown species ID for %s encountered in data file',symbol));
        end
        if verbose,
          unknown_species_ids = mytable(kinetic_data_table.(options.compound_column_name)(ind(find(cindices==0))),0)
          reaction_ids = mytable(kinetic_data_table.(options.reaction_column_name)(ind(find(cindices==0))),0)
        else
          if verbose,
            display('Set option "verbose" to see the species ID');
          end
        end
      end
      
    case 'None',   
      
  end 

  end
  
  if isempty(rindices), rindices =[]; end
  if isempty(cindices), cindices =[]; end


  %% ----------------------------------------------------------------------------
  %% Average over multiple numerical values for the same quantity

  %DISPLAY collected values
  %display(sprintf('%s - %s',symbol, related_element))
  %ind_f = find(isfinite(my_median));
  %[rindices(ind_f), cindices(ind_f),my_median(ind_f)]

  if max(rindices)>10^6, error('Weird error'); end 
  labels = rindices + 10^7 * cindices;
  labels_unique = unique(labels);
  if length(labels_unique)
    new = struct;
    for it5 = 1:length(labels_unique),
      my_ind               = find(labels_unique(it5)==labels);
      %% ALTERNATIVE: use only the last data value: my_ind=my_ind(end);
      new.rindices(it5,1)  = rindices(my_ind(1));
      new.cindices(it5,1)  = cindices(my_ind(1));
      new.my_median(it5,1) = nanmean(my_median(my_ind));
      new.my_mean(it5,1)   = nanmean(my_mean(my_ind));
      new.my_std(it5,1)    = nanmean(my_std(my_ind));
      new.my_geomstd(it5,1)    = nanmean(my_geomstd(my_ind));
      new.my_mean_ln(it5,1)= nanmean(my_mean_ln(my_ind));
      new.my_std_ln(it5,1) = nanmean(my_std_ln(my_ind));
      new.my_upper(it5,1)  = nanmin(my_upper(my_ind));
      new.my_lower(it5,1)  = nanmax(my_lower(my_ind));
    end
    ind_ok = find([new.rindices .* new.cindices] ~=0);
    rindices  = new.rindices(ind_ok);
    cindices  = new.cindices(ind_ok);
    my_median = new.my_median(ind_ok);
    my_mean   = new.my_mean(ind_ok);
    my_std    = new.my_std(ind_ok);
    my_geomstd    = new.my_geomstd(ind_ok);
    my_mean_ln= new.my_mean_ln(ind_ok);
    my_std_ln = new.my_std_ln(ind_ok);
    my_upper  = new.my_upper(ind_ok);
    my_lower  = new.my_lower(ind_ok);
  end

  %% -----------------------------------------------------------
  
  quantity_entry          = struct;
  quantity_entry.scaling  = scaling;
  quantity_entry.median = nan * zeros(ss);
  quantity_entry.mean   = nan * zeros(ss);
  quantity_entry.std    = nan * zeros(ss);
  quantity_entry.lower  = nan * zeros(ss);
  quantity_entry.upper  = nan * zeros(ss);
  switch scaling,
    case 'Multiplicative',
      quantity_entry.mean_ln = nan * zeros(ss);
      quantity_entry.std_ln  = nan * zeros(ss);
  end
  
  %% Copy values to data structure quantity_entry
  %% this doesn't work for duplicate values: only the last value is accepted
  %% also doesn't work for values that have to appear in several places in the model
  
  ind_ok = find((rindices~=0).*(cindices~=0));

  for itt = 1:length(ind_ok),
    quantity_entry.median(rindices(ind_ok(itt)),cindices(ind_ok(itt)))  = my_median(ind_ok(itt));
    quantity_entry.mean(rindices(ind_ok(itt)),cindices(ind_ok(itt)))    = my_mean(ind_ok(itt));
    quantity_entry.std(rindices(ind_ok(itt)), cindices(ind_ok(itt)))    = my_std(ind_ok(itt));
    quantity_entry.lower(rindices(ind_ok(itt)),cindices(ind_ok(itt)))   = my_lower(ind_ok(itt));
    quantity_entry.upper(rindices(ind_ok(itt)), cindices(ind_ok(itt)))  = my_upper(ind_ok(itt));
    switch scaling,
      case 'Multiplicative',
        quantity_entry.mean_ln(rindices(ind_ok(itt)), cindices(ind_ok(itt))) = my_mean_ln(ind_ok(itt));
        quantity_entry.std_ln(rindices(ind_ok(itt)), cindices(ind_ok(itt)))  = my_std_ln(ind_ok(itt));
        quantity_entry.geomstd(rindices(ind_ok(itt)), cindices(ind_ok(itt))) = my_geomstd(ind_ok(itt));
    end  
  end
  
  switch scaling,
    case 'Multiplicative',
      if flag_mean_is_given,
        [quantity_entry.mean_ln, quantity_entry.std_ln] = lognormal_normal2log(quantity_entry.mean, quantity_entry.std);
        quantity_entry.median = exp(quantity_entry.mean_ln);
      else,
        quantity_entry.mean_ln = log(quantity_entry.median);
        ii = find(isfinite(quantity_entry.std_ln));
        [quantity_entry.mean(ii), quantity_entry.std(ii)] = lognormal_log2normal(quantity_entry.mean_ln(ii), quantity_entry.std_ln(ii));
      end
      quantity_entry.geomstd = exp(quantity_entry.std_ln);
      quantity_entry.lower_ln = log(quantity_entry.lower);
      quantity_entry.upper_ln = log(quantity_entry.upper);
    case 'Additive',
      quantity_entry.geomstd = nan * quantity_entry.mean;
  end

  kinetic_data.(symbol) = quantity_entry; 
  
end

if flag_construct_empty_data_structure,
  return
end


% --------------------------------------------------------------------------------------
% make sure important fields in kinetic_data are considered

% convert Keq values from standard concentration M to standard concentration mM if necessary

standard_concentration = 'undefined';

if length(kinetic_data_sbtab),
  table_attributes = sbtab_table_get_attributes(kinetic_data_sbtab);
  if isfield(table_attributes,'StandardConcentration'),
    standard_concentration = table_attributes.StandardConcentration;
  end
end

switch standard_concentration,
  case {'M','1M','1 M'},
    ind_water = network_find_water(network);
    N_for_dissolved_compounds = network.N;
    N_for_dissolved_compounds(ind_water,:) = 0;
    molecularity_mismatch = full(sum(N_for_dissolved_compounds,1));
end

if isfield(kinetic_data,'Keq'),
  switch standard_concentration,
    case {'M','1M','1 M'},
      %% convert equilibrium constants: standard concentration M -> standard concentration mM
      kinetic_data.Keq.median   = kinetic_data.Keq.median   .* 1000.^column(molecularity_mismatch);
      kinetic_data.Keq.mean     = kinetic_data.Keq.mean     .* 1000.^column(molecularity_mismatch);
      kinetic_data.Keq.std      = kinetic_data.Keq.std      .* 1000.^column(molecularity_mismatch);
      kinetic_data.Keq.lower    = kinetic_data.Keq.lower    .* 1000.^column(molecularity_mismatch);
      kinetic_data.Keq.upper    = kinetic_data.Keq.upper    .* 1000.^column(molecularity_mismatch);
      kinetic_data.Keq.mean_ln  = kinetic_data.Keq.mean_ln   + log(1000.^column(molecularity_mismatch));
      kinetic_data.Keq.lower_ln = kinetic_data.Keq.lower_ln  + log(1000.^column(molecularity_mismatch));
      kinetic_data.Keq.upper_ln = kinetic_data.Keq.upper_ln  + log(1000.^column(molecularity_mismatch));
      if verbose,
        display('Converting equilibrium constants from a convention with standard concentration M -> standard concentration mM.');
      end
    case {'mM','1mM','1 mM'},
    case 'undefined',
      if verbose,
        display(sprintf('- (kinetic_data_load.m):\n  Thermodynamic quantities in the data file are assumed to refer to a standard concentration of mM.\n  (To change this, please add the attribute StandardConcentration=''M'' to your data file.)'));
      end
  end
end

if isfield(kinetic_data,'dmu0'),
  switch standard_concentration,
    case {'M','1M','1 M'},
    %% convert equilibrium constants: standard concentration M -> standard concentration mM
    kinetic_data.dmu0.median   = kinetic_data.dmu0.median   - RT * log(1000.^column(molecularity_mismatch));
    kinetic_data.dmu0.mean     = kinetic_data.dmu0.mean     - RT * log(1000.^column(molecularity_mismatch));
    kinetic_data.dmu0.lower    = kinetic_data.dmu0.lower    - RT * log(1000.^column(molecularity_mismatch));
    kinetic_data.dmu0.upper    = kinetic_data.dmu0.upper    - RT * log(1000.^column(molecularity_mismatch));
    if verbose,
      display('Converting standard reaction GFE from a convention with standard concentration M -> standard concentration mM.');
    end
    case {'mM','1mM','1 mM'},
    case 'undefined',
      if verbose,
        display(sprintf('- (kinetic_data_load.m):\n  Thermodynamic quantities in the data file are assumed to refer to a standard concentration of mM.\n  (To change this, please add the attribute StandardConcentration=''M'' to your data file.)'));
      end
  end
end


% --------------------------------------------------------------------------------------
% lower and upper values: make sure important fields in kinetic_data are considered

if isfield(kinetic_data,'Keq'),
  if sum(isfinite(kinetic_data.Keq.lower_ln )) == 0, 
    kinetic_data.Keq.lower_ln = log(kinetic_data.Keq.lower);
  end
  
  if sum(isfinite(kinetic_data.Keq.upper_ln )) == 0, 
    kinetic_data.Keq.upper_ln = log(kinetic_data.Keq.upper);
  end
end

if isfield(kinetic_data,'Kcatf'),
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
end


% ----------------------------------------------------------------
% if Keq values are missing, get them from dmu0 values

if isfield(kinetic_data,'Keq') * isfield(kinetic_data,'dmu0'),
  if verbose,
    display('- (kinetic_data_load.m): Completing missing Keq values based on standard Gibbs free energies of reaction (in kinetic_data_load)');
  end
  ind_missing = find(isnan(kinetic_data.Keq.median));
  kinetic_data.Keq.median(ind_missing) = exp(-1/RT * kinetic_data.dmu0.median(ind_missing));
  kinetic_data.Keq.mean_ln(ind_missing)  = -1/RT * kinetic_data.dmu0.mean(ind_missing);
  kinetic_data.Keq.std_ln(ind_missing)   =  1/RT * kinetic_data.dmu0.std(ind_missing);
  [kinetic_data.Keq.mean(ind_missing),kinetic_data.Keq.std(ind_missing)] = lognormal_log2normal(kinetic_data.Keq.mean_ln(ind_missing),kinetic_data.Keq.std_ln(ind_missing));
  kinetic_data.Keq.geomstd(ind_missing)   =  exp(kinetic_data.Keq.std_ln(ind_missing));
end
