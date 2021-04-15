function [v, v_std] = flux_data_load(network, flux_data_file, reaction_field, reaction_column, value_column, std_column)

% [v, v_std] = flux_data_load(network, flux_data_file, reaction_field, reaction_column, value_column, std_column)
%
% Load flux data (values and standard deviations for a single flux distribution) and map them to a model.
%
% o The data file must be in SBtab format
% o By default, reactions are mapped by comparing network.actions to the Reaction column in the flux file
%   (this can be changed by using the optional function arguments reaction_field and reaction_column)
% o The values must be in a column Flux, Flux:Value, or Value
% o The standard deviations must be in a column StdDev, Flux:StdDev, or they are computed from the column Flux:Quantile95

eval(default('reaction_field','''actions''', 'reaction_column','''Reaction''','value_column','[]','std_column','[]'));
  
v_data = sbtab_table_load(flux_data_file);

reaction_names = sbtab_table_get_column(v_data,reaction_column);

if isempty(value_column),
  if sbtab_table_has_column(v_data,'Flux'), value_column = 'Flux';
  elseif sbtab_table_has_column(v_data,'Flux:Value'), value_column = 'Flux:Value';
  else value_column = 'Value'; end
end

v_values = sbtab_table_get_column(v_data,value_column,1);

if isempty(std_column),
  if sbtab_table_has_column(v_data,'StdDev'), std_column = 'StdDev';
  elseif sbtab_table_has_column(v_data,'Std'), std_column = 'Std';
  elseif sbtab_table_has_column(v_data,'Flux:StdDev'), std_column = 'Flux:StdDev';
  elseif sbtab_table_has_column(v_data,'Flux:Std'), std_column = 'Flux:Std';
  else std_column = 'Value'; end
end

if length(std_column)
  v_values_std = sbtab_table_get_column(v_data,std_column,1);
else
  v_values_std = [];
  if length(sbtab_table_get_column(v_data,'Flux:Quantile95',1,0)),
    %% 95% of all values are within 1.96 sigma)
    v_values_std   = 1/[2*1.96] * sbtab_table_get_column(v_data,'Flux:Quantile95',1); % 
  end
end

% --------------------------------------------------

v           = nan * ones(size(network.(reaction_field)));
ll          = label_names(network.(reaction_field),reaction_names);
v(find(ll)) = v_values(ll(find(ll)));
v(strcmp('',network.(reaction_field))) = nan; 

if sum(isnan(v)), 
  warning('Some flux data are missing'); 
end

if length(v_values_std),
  v_std  = nan * ones(size(network.(reaction_field)));
  v_std(find(ll))= v_values_std(ll(find(ll)));
  v_std(strcmp('',network.(reaction_field))) = nan; 
else 
  v_std = [];
end
