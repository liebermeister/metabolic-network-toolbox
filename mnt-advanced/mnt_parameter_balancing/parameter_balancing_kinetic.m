function [r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std] = parameter_balancing_kinetic(network, kinetic_data, options);

%  [r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std] = parameter_balancing_kinetic(network, kinetic_data, options);
%
% Determine consistent kinetic parameter set by parameter balancing
%
% This is a convenience function for parameter balancing, with settings adjusted for kinetic constants
%
% Input
%   network         Metabolic network data structure
%   kinetic_data    Kinetic_data used
%   options         Parameter balancing options used
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
% options: struct with options; for default values, see 'parameter_balancing_default_options'
% this function uses (potentially) the following options
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
%
% The code assumes that the network structure contains KEGG IDs and uses SBtab data files 
% These files can contain data about the following biochemical quantities: 
%   standard chemical potential
%   Michaelis constant
%   activation constant  
%   inhibitory constant
%   equilibrium constant 
%   substrate catalytic rate constant 
%   product catalytic rate constant
%
% If "kinetic_data" is a string, it is taken to be a filename; otherwise is is taken to be a data structure
%
% The standard reaction directions in the model have to follow the convention in KEGG!!

  
% ------------------------------------------------------------------------
% Initialise
  
eval(default('kinetic_data', '[]', 'options', 'struct'));

options = join_struct(parameter_balancing_default_options,options);

[nm,nr] = size(network.N);


% ----------------------------------------------------------------
% Lists of quantities to be considered
% (consider using 'parameter_balancing_quantities' instead)

basic_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant'}';

model_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant', 'Michaelis constant product'}';

data_quantities   = {'standard Gibbs energy of reaction', 'standard chemical potential','Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'}';


% ------------------------------------------------------------------------
% Load and modify parameter priors

parameter_prior = parameter_balancing_prior([],options.parameter_prior_file);

parameter_prior = pb_parameter_prior_adjust(parameter_prior, options); 


% ----------------------------------------------------------------
% Load and modify kinetic data

if isstr(kinetic_data),
  %% Load data from file
  kinetic_data = data_integration_load_kinetic_data(data_quantities, [], network, kinetic_data, 0, 1);
elseif isempty(kinetic_data),  
  %% If necessary, create empty kinetic_data structure
  kinetic_data = data_integration_load_kinetic_data(data_quantities, [], network, [], 0, 1, options.reaction_column_name, options.compound_column_name);
end  

kinetic_data_orig = kinetic_data;

kinetic_data      = pb_kinetic_data_adjust(kinetic_data, parameter_prior, network, options);


% ----------------------------------------------------------------
% Run parameter balancing

network.kinetics                   = set_kinetics(network, 'cs');
task                               = parameter_balancing_task(network, kinetic_data, parameter_prior, model_quantities, basic_quantities);
res                                = parameter_balancing_calculation(task, parameter_prior, options);
[r,r_mean,r_std,r_orig,r_samples]  = parameter_balancing_output(res, kinetic_data_orig, options);


