function options = parameter_balancing_default_options()
  
% Options and their default values:
%
%   options.kinetics             = 'cs';
%   options.parametrisation      = 'catalytic rate constant'; 
%   options.reaction_units       = 'concentration per time'; % 'amount per time'
%   options.include_metabolic    = 1; % consider concentrations as quantities in parameter balancing
%   options.enzymes_explicit     = 1; % do not consider maximal velocities as quantities in parameter balancing
%   options.kcat_prior_median    = [];
%   options.kcat_prior_log10_std = []
%   options.KM_lower             = [];
%   options.Keq_upper            = [];
%   options.kcat_lower           = [];
%   options.kcatr_lower          = [];
%   options.kcat_upper           = []; 
%   options.GFE_fixed            = 0
%   options.parameter_prior_file = []; % => use default file, see biochemical_parameter_prior
%   options.use_pseudo_values    = 1;
%   options.n_samples            = 0;
%   options.Keq_given            = [];
%   options.enforce_flux_directions = 0;
%   options.postprocessing_enforce_fluxes = 0;

options.kinetics             = 'cs';
options.parametrisation      = 'catalytic rate constant'; 
options.reaction_units       = 'concentration per time'; % 'amount per time'
options.enzymes_explicit     = 1;
options.include_metabolic    = 1;
options.kcat_prior_median    = [];
options.kcat_prior_log10_std = []; 
options.kcat_usage           = 'use';
options.KM_lower             = [];
options.Keq_upper            = [];
options.kcat_lower           = []; 
options.kcatr_lower          = []; 
options.kcat_upper           = []; 
options.GFE_fixed            = []; 
options.parameter_prior_file = [];
options.use_pseudo_values    = 0;
options.n_samples            = 0; 
options.Keq_given            = []; 
options.postprocessing_adjustment_to_fluxes = 0;
