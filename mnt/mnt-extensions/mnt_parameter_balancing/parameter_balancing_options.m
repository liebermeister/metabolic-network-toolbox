function options = parameter_balancing_options()

% PARAMETER_BALANCING_OPTIONS Options for parameter balancing - Data structure with default values
%
% options = parameter_balancing_options()
%
% General options
%
%   options.parameter_prior_file         = [];                        default prior file (see parameter_balancing_prior.m)
%   options.kinetics                     = 'cs';                      kinetics for SBML output model (default common modular rate law)
%   options.parametrisation              = 'catalytic rate constant'; type of model parameterisation
%   options.reaction_units               = 'concentration per time';  physical unit. alternative: 'amount per time'
%   options.GFE_fixed                    = 0                          fixed GFE values          
%   options.insert_pseudo_values         = 0;                         replace missing data values by pseudo values
%   options.use_bounds_from_prior_table  = 1;                         flag: use bounds
%   options.include_metabolic            = 1;                         flag: consider concentrations in parameter balancing
%   options.enzymes_explicit             = 1;                         flag: disregard maximal velocities in parameter balancing
%   options.enforce_flux_directions      = 0;                         flag: set constraints to realise predefined flux directions
%   options.adjust_to_fluxes             = 0;                         flag: apply postprocessing to yield predefined fluxes (see parameter_balancing_output)
%   options.v                            = [];                        flux vector used by options "enforce_flux_directions" and "adjust_to_fluxes"
%   options.ignore_all_constraints       = 0;                         flag: run parameter balancing without constraints
%   options.n_samples                    = 0;                         0: no sampling; n>1: generate n sampled parameter sets
%   options.preferred_data_element_ids   = 'id';                      data columns to be used for mapping (default: "Reaction" or "Compound") % alternatives 'sbml', 'kegg'
%   options.use_sbml_ids                 = 0;                           (redundant with preferred_data_element_ids; field used internally)
%   options.use_kegg_ids                 = 0;                           (redundant preferred_data_element_ids; field used internally)


  
  
% Options to override bounds (given in prior or data file)
%
%   options.conc_min_default                  = 0.001;                concentration default lower bound 0.001 mM
%   options.conc_max_default                  = 10;                   concentration default upper bound 10 mM
%   options.conc_min                          = 0.00001;              concentration lower bound default (mM)
%   options.conc_max                          = 100;                  concentration upper bound default (mM)
%   options.kcat_prior_median                 = [];                   predefined kcat prior median value
%   options.kcat_prior_log10_std              = []                    predefined kcat prior log 10 std dev
%   options.kcat_usage                        = 'use';                flag: use given kcat values 
%   options.KM_lower                          = []; % mM              KM value default lower bound
%   options.KM_upper                          = []; % mM              KM value default upper bound 
%   options.Keq_upper                         = [];                   Keq value default lower bound
%   options.Keq_given                         = [];                   Keq values (given vector)
%   options.kcat_lower                        = []; % 1/s             kcat value default lower bound
%   options.kcatr_lower                       = []; % 1/s             reverse kcat value default lower bound
%   options.kcat_upper                        = []; % 1/s             kcat value default upper bound
                                             
% Optionns for wrapper functions             
%
%   options.use_pseudo_values                 = 1;                    flag: use pseudo values
%   options.use_data                          = 1;                    use data table
%   options.flag_check                        = 1;                    check results and show diagnostic plots
%
% Extra options used in 'parameter_balancing_thermodynamic'

%   options.ind_water                         = []                    indices of metabolites representing water
%   options.set_water_conc_to_one             = 1;                    flag: preprocessing: set water concentration to one
%   options.data_refer_to_molar               = 0;                    flag: concentrations in molar (default: millimolar)
%   options.A_max                             = 1000;                 A lower bound (default) (kJ/mol)
%   options.A_min                             = 0.5;                  A lower bound (default) (kJ/mol)
%   options.A_mean                            = nan;                  A mean (default) (kJ/mol)
%   options.A_std                             = std;                  A mean (default) (kJ/mol)
%   options.A_lower                           = nan * ones(nr,1);     A lower bounds (vector)
%   options.A_upper                           = nan * ones(nr,1);     A upper bounds (vector)
%   options.c_min                             = nan                   (vector)
%   options.c_max                             = nan                   (vector)
%   options.c_fix                             = nan                   (vector)
%   options.c_fix_strict                      = 0;                    strictly fix concentrations (otherwise, allow for a range)!
%   options.u_max                             = 1000;                 (mM)
%   options.u_min                             = 0.01;                 (mM)
%   options.variability                       = 2;                    variability of known concentrations
%   options.sigma_mu0                         = 3;                    error of mu0 values (kJ/mol); 3 for alberty data
%   options_default.parameter_prior_filename  = [];                   file containing the prior table to be used
%   options_default.test_v_for_being_feasible = 1;                    run previous test for feasible v
%   options_default.fix_Keq_in_sampling       = 0;                    flag
%
% CURRENTLY NOT USED (used in python code only);
%   options.ph	                
%   options.temperature	        
%   options.overwrite_kinetics	
%   options.cell_volume	        
%   options.enzyme_prefactor	
%   options.default_inhibition	
%   options.default_activation	
%   options.model_name	        
%   options.boundary_value
  
options.parameter_prior_file        = [];
options.kinetics                    = 'cs'; % common modular rate law
options.parametrisation             = 'catalytic rate constant'; 
options.reaction_units              = 'concentration per time'; % alternative 'amount per time'
options.conc_min_default            = 0.001; % mM
options.conc_max_default            = 10;    % mM
options.conc_min                    = 0.00001; %(mM)
options.conc_max                    = 100;     %(mM)
options.kcat_prior_median           = [];
options.kcat_prior_log10_std        = []; 
options.KM_lower                    = [];
options.KM_upper                    = [];
options.Keq_upper                   = [];
options.Keq_given                   = []; 
options.kcat_lower                  = []; 
options.kcatr_lower                 = []; 
options.kcat_upper                  = []; 
options.GFE_fixed                   = []; 
options.kcat_usage                  = 'use';
options.use_pseudo_values           = 1;
options.insert_pseudo_values        = 0;
options.use_bounds_from_prior_table = 1;
options.include_metabolic           = 1;
options.enzymes_explicit            = 1;
options.enforce_flux_directions     = 0;
options.adjust_to_fluxes            = 0;
options.ignore_all_constraints      = 0;
options.n_samples                   = 0; 
options.preferred_data_element_ids  = 'id'; % alternatives 'sbml_id', 'kegg_id'
options.use_sbml_ids                = 0;
options.use_kegg_ids                = 0;

options.flag_check                  = 1;
options.use_data                    = 1;

options.ind_water                   = [];
options.set_water_conc_to_one       = 1;
options.data_refer_to_molar         = 0;
options.sigma_mu0                   = 3;%    error of mu0 values (kJ/mol); 3 for alberty data
options.A_max                       = 100;
options.A_min                       = 0.5;  
options.c_min                       = [];
options.c_max                       = [];
options.c_fix                       = [];
options.c_fix_strict                = 0;
options.A_fix                       = [];
options.A_lower                     = [];
options.A_upper                     = [];
options.virtual_reactions           = {};
options.parameter_prior_filename    = [];
options.test_v_for_being_feasible   = 1;
options.fix_Keq_in_sampling         = 0;

