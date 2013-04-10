function [task, prior] = parameter_balancing_task(network, kinetic_data, quantity_info, model_quantities, basic_quantities)

% Prepare vectors and matrices for parameter balancing
%
% Mandatory inputs: network, kinetic_data
%
% All q and x values are given in natural scaling!
%
% 'prior' has the same data format as kinetic_data (for export to table files)

eval(default('quantity_info','[]'));

if isempty(quantity_info),
  quantity_info =   data_integration_load_quantity_info;
end

[nr,nm,nx,ind_KM,ind_KA,ind_KI,nKM,nKA,nKI] = network_numbers(network);

if ~exist('model_quantities','var'),
  [model_quantities, basic_quantities] = parameter_balancing_quantities(quantity_info,network); 
end


% -------------------------------------------------------------
% initialise data structure 'task'

% for which quantities are data available?

data_quantities = quantity_info.QuantityType(label_names(fieldnames(kinetic_data),quantity_info.Symbol));


task.network          = network;
task.model_quantities = model_quantities;
task.basic_quantities = basic_quantities;
task.data_quantities  = data_quantities;
task.q.names          = {};
task.q.scaling        = {};
task.q.indices        = struct;
task.q.prior.mean     = [];
task.q.prior.std      = [];


% -------------------------------------------------------------
% build description of basic quantities

num_basic = parameter_balancing_quantity_numbers(basic_quantities,quantity_info,network);

for it = 1:length(basic_quantities),
  my_quantity = basic_quantities{it};
  ind = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_scaling      = quantity_info.Scaling{ind};
  my_symbol       = quantity_info.Symbol{ind};
  my_prior_median = cell_string2num(quantity_info.PriorMedian(ind));
  my_prior_std    = cell_string2num(quantity_info.PriorStd(ind));
  my_rel_element  = quantity_info.RelatedElement{ind};
  
  %% make data structure "prior"

  prior.(my_symbol).scaling = my_scaling;
  prior.(my_symbol).median  = my_prior_median;

  switch my_rel_element, 
    case 'Species',   ee = ones(length(network.metabolites), 1);
    case 'Reaction',  ee = ones(length(network.actions), 1);
    case 'Reaction/Species', ee = ones(length(network.actions), length(network.metabolites));
  end
  
  if strcmp(my_scaling, 'Logarithmic'),
    prior.(my_symbol).mean_ln = log(my_prior_median) * ee;
    prior.(my_symbol).std_ln  = log(10) * my_prior_std * ee;
    [prior.(my_symbol).mean, prior.(my_symbol).std] = lognormal_log2normal(prior.(my_symbol).mean_ln, prior.(my_symbol).std_ln);
  else,
    prior.(my_symbol).mean = my_prior_median * ee;
    prior.(my_symbol).std  = my_prior_std * ee;
  end
  
  
  if strcmp(my_scaling, 'Logarithmic'),
      my_prior_median = log(my_prior_median);
      my_prior_std    = log(10) * my_prior_std;
  end
  my_indices = length(task.q.prior.mean) + [1:num_basic(it)]';
  task.q.names(my_indices,:)   = repmat(basic_quantities(it), num_basic(it),1);
  task.q.scaling(my_indices,:) = repmat({my_scaling}, num_basic(it),1);
  task.q.indices.(my_symbol)   = my_indices;
  task.q.prior.mean            = [task.q.prior.mean; my_prior_median * ones(num_basic(it),1)];
  task.q.prior.std             = [task.q.prior.std ; my_prior_std    * ones(num_basic(it),1)];

end

task.q.numbers = num_basic;


% -------------------------------------------------------------
% build description of model quantities

task.xmodel.names    = {};
task.xmodel.scaling  = {};
task.xmodel.indices  = struct;

[nr,nm,nx,ind_KM,ind_KA,ind_KI,nKM,nKA,nKI] = network_numbers(network);

num_model = parameter_balancing_quantity_numbers(model_quantities,quantity_info,network);

for it = 1:length(model_quantities),
  my_quantity     = model_quantities{it};
  ind             = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_scaling      = quantity_info.Scaling{ind};
  my_symbol       = quantity_info.Symbol{ind};
  my_indices      = length(task.xmodel.names) + [1:num_model(it)]';
  task.xmodel.names(my_indices,:)   = repmat(model_quantities(it), num_model(it),1);
  task.xmodel.scaling(my_indices,:) = repmat({my_scaling}, num_model(it),1);
  task.xmodel.indices.(my_symbol)   = my_indices;      
end

task.xmodel.numbers = num_model;

Q_model = parameter_balancing_construct_Q_matrix(model_quantities, basic_quantities, quantity_info, network);

task.Q_xmodel_q = Q_model;


% -------------------------------------------------------------
% build description of all potential relevant quantities
% (for pseudo values and general constraints)

all_quantities = unique([model_quantities; basic_quantities; data_quantities]);
num_all        = parameter_balancing_quantity_numbers(all_quantities,quantity_info,network);

x_all.names       = {};
x_all.scaling     = {};
x_all.pseudo      = struct;
x_all.pseudo.mean = [];
x_all.pseudo.std  = [];
x_all.upper_nat   = [];
x_all.lower_nat   = [];
x_all.indices     = struct;

for it = 1:length(all_quantities),

  my_quantity                 = all_quantities{it};
  ind                         = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_scaling                  = quantity_info.Scaling{ind};
  my_symbol                   = quantity_info.Symbol{ind};
  my_indices                  = length(x_all.names) + [1:num_all(it)]';
  x_all.names(my_indices,:)   = repmat(all_quantities(it), num_all(it),1);
  x_all.scaling(my_indices,:) = repmat({my_scaling}, num_all(it),1);
  x_all.indices.(my_symbol)   = my_indices;

  my_all_upper    = eval(quantity_info.UpperBound{ind});
  my_all_lower    = eval(quantity_info.LowerBound{ind});
  my_pseudo_median= cell_string2num(quantity_info.PriorMedian(ind));
  my_pseudo_std   = cell_string2num(quantity_info.PriorStd(ind));

  switch my_scaling, case 'Logarithmic',
    my_all_lower = log(my_all_lower);
    my_all_upper = log(my_all_upper);
    my_pseudo_median = log(my_pseudo_median);
    my_pseudo_std    = log(10) * my_pseudo_std;
  end
  
  x_all.pseudo.mean = [x_all.pseudo.mean; my_pseudo_median * ones(num_all(it),1)];
  x_all.pseudo.std  = [x_all.pseudo.std ; my_pseudo_std    * ones(num_all(it),1)];
  x_all.upper_nat   = [x_all.upper_nat;   my_all_upper * ones(num_all(it),1)];
  x_all.lower_nat   = [x_all.lower_nat;   my_all_lower * ones(num_all(it),1)];

end


% -------------------------------------------------------------
% build description of data quantities

x_data_mean  = [];
x_data_std   = [];
x_data_lower = [];
x_data_upper = [];

x_data.names   = {};
x_data.scaling = {};
num_data     = parameter_balancing_quantity_numbers(data_quantities,quantity_info,network);

for it = 1:length(data_quantities),
  
  my_quantity                   = data_quantities{it};
  ind                           = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_scaling                    = quantity_info.Scaling{ind};
  my_symbol                     = quantity_info.Symbol{ind};
  my_indices                    = length(x_data.names) + [1:num_data(it)]';
  x_data.names(my_indices,:)    = repmat(data_quantities(it), num_data(it),1);
  x_data.scaling(my_indices,:)  = repmat({my_scaling}, num_data(it),1);
  x_data.indices.(my_symbol)    = my_indices;

  switch my_scaling,
    
    case 'Original',
      my_x_mean  = kinetic_data.(my_symbol).mean;
      my_x_std   = kinetic_data.(my_symbol).std;
      my_x_lower = kinetic_data.(my_symbol).lower;
      my_x_upper = kinetic_data.(my_symbol).upper;
    
    case 'Logarithmic',
      my_x_mean  = kinetic_data.(my_symbol).mean_ln;
      my_x_std   = kinetic_data.(my_symbol).std_ln;
      my_x_lower = kinetic_data.(my_symbol).lower_ln;
      my_x_upper = kinetic_data.(my_symbol).upper_ln;
  end
  
  switch my_symbol, 
    
    case 'KM', 
      my_x_mean  = my_x_mean(ind_KM);
      my_x_std   = my_x_std(ind_KM);
      my_x_lower = my_x_lower(ind_KM);
      my_x_upper = my_x_upper(ind_KM);
    
    case 'KA', 
      my_x_mean  = my_x_mean(ind_KA);
      my_x_std   = my_x_std(ind_KA);
      my_x_lower = my_x_lower(ind_KA);
      my_x_upper = my_x_upper(ind_KA);
    
    case 'KI', 
      my_x_mean  = my_x_mean(ind_KI);
      my_x_std   = my_x_std(ind_KI);
      my_x_lower = my_x_lower(ind_KI);
      my_x_upper = my_x_upper(ind_KI);
  
  end

  x_data_mean  = [x_data_mean;  column(my_x_mean(:))  ];
  x_data_std   = [x_data_std ;  column(my_x_std(:))   ];
  x_data_lower = [x_data_lower; column(my_x_lower(:)) ];
  x_data_upper = [x_data_upper; column(my_x_upper(:)) ];

end


% -------------------------------------------------------------------
% construct the matrix for all relevant quantities

Q_all  = parameter_balancing_construct_Q_matrix(all_quantities, basic_quantities, quantity_info, network);
Q_data = parameter_balancing_construct_Q_matrix(data_quantities, basic_quantities, quantity_info, network);


% -------------------------------------------------------------------
% reduce vectors and matrix to existing data

ind_rel = find(isfinite(x_data_mean));
task.Q_xdata_q     = Q_data(ind_rel,:);

task.xdata.names   = x_data.names(ind_rel);
task.xdata.scaling = x_data.scaling(ind_rel);
task.xdata.mean    = x_data_mean(ind_rel);
task.xdata.std     = x_data_std(ind_rel);
task.xdata.mean    = x_data_mean(ind_rel);
task.xdata.indices = struct;
task.xdata.numbers = [];

for it = 1:length(data_quantities),
  my_quantity = data_quantities{it};
  ind         = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_symbol   = quantity_info.Symbol{ind};
  my_ind      = find(strcmp(my_quantity,task.xdata.names));
  if length(my_ind), 
    task.xdata.indices.(my_symbol) = my_ind;
    task.xdata.numbers = [task.xdata.numbers; length(my_ind)];
  end
end


task.xall       = x_all;
task.Q_xall_q   = Q_all;



% -------------------------------------------------------------------
% reduce vectors and matrix to existing data, for lower bounds

ind_rel = isfinite(x_data_lower);

task.xlower.names     = x_data.names(ind_rel);
task.xlower.scaling   = x_data.scaling(ind_rel);
task.xlower.value_nat = x_data_lower(ind_rel);
task.xlower.indices   = struct;
task.xlower.numbers   = [];

for it = 1:length(data_quantities),
  my_quantity = data_quantities{it};
  ind         = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_symbol   = quantity_info.Symbol{ind};
  my_ind      = find(strcmp(my_quantity,task.xlower.names));
  if length(my_ind), 
    task.xlower.indices.(my_symbol) = my_ind;
    task.xlower.numbers = [task.xlower.numbers; length(my_ind)];
  end
end

task.Q_xlower_q   = Q_data(ind_rel,:);


% -------------------------------------------------------------------
% reduce vectors and matrix to existing data, for upper bounds

ind_rel = isfinite(x_data_upper);

task.xupper.names   = x_data.names(ind_rel);
task.xupper.scaling = x_data.scaling(ind_rel);
task.xupper.value_nat = x_data_upper(ind_rel);
task.xupper.indices = struct;
task.xupper.numbers = [];

for it = 1:length(data_quantities),
  my_quantity = data_quantities{it};
  ind         = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_symbol   = quantity_info.Symbol{ind};
  my_ind      = find(strcmp(my_quantity,task.xupper.names));
  if length(my_ind),
    task.xupper.indices.(my_symbol) = my_ind;
    task.xupper.numbers = [task.xupper.numbers; length(my_ind)];
  end
end

task.Q_xupper_q = Q_data(ind_rel,:);

% --------------------------------------------------------------------

if find(isnan(task.xdata.std)), 
  warning('Some standard deviations are missing; replacing them by 1.');
  task.xdata.std(find(isnan(task.xdata.std))) = 1; 
end
