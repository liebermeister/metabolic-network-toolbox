function [model_quantities, basic_quantities, data_quantities] = parameter_balancing_quantities(quantity_info,network,options)

% [model_quantities, basic_quantities, data_quantities] = parameter_balancing_relevant_quantities(quantity_info,network,options)
% 
% Lists of 
%  - quantities relevant for a model (model_quantities, list of strings),
%  - relevant basic quantities in a model (basic_quantities)
%  - possible quantities in kinetic data that might be relevant to estimate the basic quantities
%
% Options
%  options.include_metabolic (default 0) 
%  options.enzymes_explicit  (default 1)
%  options.parametrisation   (default 'catalytic rate constant');

eval(default('quantity_info','[]','options','struct'));

if isempty(quantity_info),
  quantity_info =   data_integration_load_quantity_info;
end

options_default = struct('include_metabolic',0,'enzymes_explicit',1,'parametrisation','catalytic rate constant');
options = join_struct(options_default,options);

% ------------------------------------------------
% restrict quantity table to quantities relevant for this model

% which quantities are needed in the model ? 

quantities = quantity_info.QuantityType; 

keep = ones(size(quantities));

if options.include_metabolic ==0,
  keep = keep .* double(~strcmp(quantity_info.Constant,'Dynamic'));
end

if options.enzymes_explicit,
  keep(find(strcmp(quantity_info.QuantityType,'forward maximal velocity')))   = 0;
  keep(find(strcmp(quantity_info.QuantityType,'reverse maximal velocity')))   = 0;
end

keep_model = keep;
keep_data  = keep;
keep_basic = keep;

keep_model(find(strcmp(quantity_info.QuantityType,'forward enzyme mass action term')))   = 0;
keep_model(find(strcmp(quantity_info.QuantityType,'reverse enzyme mass action term'))) = 0;

switch options.parametrisation,
  case 'equilibrium constant',
    keep_model(find(strcmp(quantity_info.QuantityType,'standard chemical potential')))   = 0;
    keep_model(find(strcmp(quantity_info.QuantityType,'substrate catalytic rate constant'))) = 0;
    keep_model(find(strcmp(quantity_info.QuantityType,'product catalytic rate constant'))) = 0;
  case 'standard chemical potential',
    keep_model(find(strcmp(quantity_info.QuantityType,'equilibrium constant')))   = 0;
    keep_model(find(strcmp(quantity_info.QuantityType,'substrate catalytic rate constant'))) = 0;
    keep_model(find(strcmp(quantity_info.QuantityType,'product catalytic rate constant'))) = 0;
  case 'catalytic rate constant',
    %% keep eq. const in order to impose pseudo values 
    %%keep_model(find(strcmp(quantity_info.QuantityType,'equilibrium constant')))   = 0;
    keep_model(find(strcmp(quantity_info.QuantityType,'standard chemical potential')))   = 0;
    keep_model(find(strcmp(quantity_info.QuantityType,'catalytic rate constant geometric mean')))    = 0;
end

keep_data(find(strcmp(quantity_info.QuantityType,'catalytic rate constant geometric mean')))   = 0;

keep_basic = keep_basic .* double(strcmp(quantity_info.Dependence,'Basic'));

% basic parameters could be further reduced to those on which any of the 
% model or data parameters depend ...

data_quantities  = quantities(find(keep_data)); 
model_quantities = quantities(find(keep_model)); 
basic_quantities = quantities(find(keep_basic));
