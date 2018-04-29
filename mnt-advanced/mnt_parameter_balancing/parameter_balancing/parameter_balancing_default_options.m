function options = parameter_balancing_default_options()

% PARAMETER_BALANCING_DEFAULT_OPTIONS Options for parameter balancing - empty data structure
%
% options = parameter_balancing_default_options()
%
%   options.parameter_prior_file = []; % => use default file, see parameter_balancing_prior.m
%   options.kinetics             = 'cs';
%   options.parametrisation      = 'catalytic rate constant'; 
%   options.reaction_units       = 'concentration per time'; % 'amount per time'
%   options.conc_min_default     = 0.001; % mM
%   options.conc_max_default     = 10;  % mM
%   options.conc_min             = [];
%   options.conc_max             = [];
%   options.kcat_prior_median    = [];
%   options.kcat_prior_log10_std = []
%   options.kcat_usage           = 'use';
%   options.KM_lower             = []; % mM
%   options.KM_upper             = []; % mM
%   options.Keq_upper            = [];
%   options.Keq_given            = [];
%   options.kcat_lower           = []; % 1/s
%   options.kcatr_lower          = []; % 1/s
%   options.kcat_upper           = []; % 1/s
%   options.GFE_fixed            = 0
%   options.use_pseudo_values    = 1;
%   options.insert_pseudo_values = 0; % replace missing data values by pseudo values
%   options.use_bounds_from_quantity_table = 1;
%   options.include_metabolic    = 1; % consider concentrations as quantities in parameter balancing
%   options.enzymes_explicit     = 1; % do not consider maximal velocities as quantities in parameter balancing
%   options.enforce_flux_directions = 0;
%   options.postprocessing_adjustment_to_fluxes = 0;
%   options.n_samples            = 0;
%
% Currently not used (used in python code only);
%   options.use_pseudos	True
%   options.ph	7
%   options.temperature	300
%   options.overwrite_kinetics	True
%   options.cell_volume	43
%   options.enzyme_prefactor	True
%   options.default_inhibition	complete
%   options.default_activation	complete
%   options.model_name	outputname
%   options.boundary_values	ignore
  
options.parameter_prior_file = [];
options.kinetics             = 'cs';
options.parametrisation      = 'catalytic rate constant'; 
options.reaction_units       = 'concentration per time'; % 'amount per time'
options.conc_min_default     = 0.001; % mM
options.conc_max_default     = 10;    % mM
options.conc_min             = [];
options.conc_max             = [];
options.kcat_prior_median    = [];
options.kcat_prior_log10_std = []; 
options.kcat_usage           = 'use';
options.KM_lower             = [];
options.KM_upper             = [];
options.Keq_upper            = [];
options.Keq_given            = []; 
options.kcat_lower           = []; 
options.kcatr_lower          = []; 
options.kcat_upper           = []; 
options.GFE_fixed            = []; 
options.use_pseudo_values    = 1;
options.insert_pseudo_values = 0;
options.use_bounds_from_quantity_table = 1;
options.include_metabolic    = 1;
options.enzymes_explicit     = 1;
options.enforce_flux_directions = 0;
options.postprocessing_adjustment_to_fluxes = 0;
options.n_samples            = 0; 
