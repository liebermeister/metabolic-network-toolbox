function [r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std,r_geom_mean,r_geom_std] = parameter_balancing_kinetic(network, kinetic_data, options);

% PARAMETER_BALANCING_KINETIC Determine a consistent kinetic parameter set
%
% [r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std,r_geom_mean,r_geom_std] = parameter_balancing_kinetic(network, kinetic_data, options);
%
% Wrapper function for parameter balancing, with matlab data structures (for model and data) as input arguments
%
% Input
%   network         Metabolic network data structure
%   kinetic_data    Kinetic_data used
%   options         Parameter balancing options used (see 'parameter_balancing_options')
%
% Output
%   r               Kinetic constants (posterior mode values, respecting the linear constraints)
%   r_orig          Original kinetic constants (used as input in parameter balancing)
%   r_mean          Kinetic constants (posterior means, ignoring the  linear constraints)
%   r_std           Kinetic constants (posterior standard deviations, ignoring the  linear constraints)
%   r_samples       Kinetic constants sampled from the posterior
%   kinetic_data    Kinetic_data used
%   parameter_prior Prior distributions used
% 
% This function uses (potentially) the options
%
%   options.flag_given_kinetics
%   options.reaction_column_name  (only if no kinetic data are given)
%   options.compound_column_name  (only if no kinetic data are given)
%   options.kcat_usage            {'use','none','forward'} (default: 'use')
%   options.kcat_prior_median
%   options.kcat_prior_log10_std
%   options.kcat_lower
%   options.kcatr_lower
%   options.kcat_upper
%   options.KM_lower
%   options.Keq_upper
%   options.parameter_prior_file
%   options.GFE_fixed
%   options.use_pseudo_values
%   options.fix_Keq_in_sampling
%   options.adjust_to_fluxes
%   options.v
%
% The function assumes that the network structure contains KEGG IDs
% The standard reaction directions in the model must follow the convention in KEGG
%
% If "kinetic_data" is a string, it is used as an SBtab data filename; otherwise is is used as a data structure
%
% The SBtab data file can contain values for the following biochemical quantities: 
%   - standard chemical potential
%   - Michaelis constant
%   - activation constant  
%   - inhibitory constant
%   - equilibrium constant 
%   - substrate catalytic rate constant 
%   - product catalytic rate constant
%
% If the option "adjust_to_fluxes" is set, fluxes must be given in the options field "v"

  
% ------------------------------------------------------------------------
% Initialise some variables

eval(default('kinetic_data', '[]', 'options', 'struct'));

options = parameter_balancing_update_options(join_struct(parameter_balancing_options,options));

[nm,nr] = size(network.N);


% ----------------------------------------------------------------
% Lists of quantities to be considered 
% (function 'parameter_balancing_quantities' can be used instead)

basic_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant'}';

pseudo_quantities  = {'equilibrium constant','reaction affinity'}';

model_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant', 'Michaelis constant product'}';

data_quantities   = {'standard Gibbs energy of reaction', 'standard chemical potential','Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'}';

if options.include_metabolic,
  basic_quantities  = [ basic_quantities;  {'concentration','concentration of enzyme'}' ];
  pseudo_quantities = [ pseudo_quantities; ];
  model_quantities  = [ model_quantities;  {'Michaelis constant product', 'concentration','reaction affinity','concentration of enzyme'}' ];
  data_quantities   = [ data_quantities;   {'concentration','reaction affinity','concentration of enzyme'}' ];
end

% ------------------------------------------------------------------------
% Load and modify parameter priors

parameter_prior = parameter_balancing_prior([],options.parameter_prior_file);
parameter_prior = pb_parameter_prior_adjust(parameter_prior, options); 


% ----------------------------------------------------------------
% Load and preprocess kinetic data
%   if kinetic_data is empty   : Create empty_data structure
%   if kinetic_data is a string: Load data

if isempty(kinetic_data),  
  kinetic_data = data_integration_load_kinetic_data(data_quantities, [], network, [], 1, 0, options.reaction_column_name, options.compound_column_name);
elseif isstr(kinetic_data),
  kinetic_data = data_integration_load_kinetic_data(data_quantities, [], network, kinetic_data, options.use_sbml_ids, options.use_kegg_ids);
end

kinetic_data_orig = kinetic_data;

kinetic_data = pb_kinetic_data_adjust(kinetic_data, parameter_prior, network, options);

% display the adjusted data
% data_integration_display_kinetic_data(kinetic_data,network);

% ----------------------------------------------------------------
% Run parameter balancing

network.kinetics  = set_kinetics(network, 'cs');

task   = parameter_balancing_task(network, kinetic_data, parameter_prior, model_quantities, basic_quantities, pseudo_quantities);
result = parameter_balancing_calculation(task, parameter_prior, options);

[r,r_mean,r_std,r_geom_mean,r_geom_std,r_orig,r_samples] = parameter_balancing_output(result, kinetic_data_orig, options);
