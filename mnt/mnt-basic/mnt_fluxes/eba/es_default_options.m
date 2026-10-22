function [es_options, es_constraints] = es_default_options(N)

% ES_DEFAULT_OPTIONS - Set default values in 'es_options' and 'es_constraints'
%
% [es_options, es_constraints] = es_default_options(N)
% 
% Initialise structures 'es_options' and 'es_constraints' used to define tasks in elasticity sampling
% 
% Function arguments: 
% N can be the (complete) stoichiometric matrix or a 'network' data structure from the Metabolic Network Toolbox)
% 
% Function outputs:
% 
% Structure 'es_constraints': directives for defining the constraints in model sampling
%
% Field name                        Size    Meaning
% ----------                        ----    -------
% es_constraints.v_fix              nr x 1  fluxes predefined values (overrides vmin and vmax)
% es_constraints.v_min              nr x 1  fluxes lower bounds (overrides -vmax)
% es_constraints.v_max              nr x 1  fluxes upper bounds
% es_constraints.v_sign             nr x 1  fluxes signs (overrides vmin and vmax)
% es_constraints.v_mean             nr x 1  fluxes data (can also be a matrix)
% es_constraints.v_std              nr x 1  fluxes data error bars (can also be a matrix)
% es_constraints.Keq_fix            nr x 1  equilibrium constants predefined values 
% es_constraints.ext_sign           nx x 1  signs of production rates of external metabolites
% es_constraints.log_c_fix          nm x 1  ln metabolite levels predefined values
% es_constraints.log_c_min          nm x 1  ln metabolite levels lower bounds 
% es_constraints.log_c_max          nm x 1  ln metabolite levels upper bounds 
% es_constraints.log_c_mean         nm x 1  ln metabolite levels data
% es_constraints.log_c0_mean        nm x 1  ln equilibrium metabolite levels data
% es_constraints.log_c0_std         nm x 1  ln equilibrium metabolite levels data error bars
% es_constraints.mu0_fix            nm x 1  std chemical potentials predefined values 
% es_constraints.mu0_min            nm x 1  std chemical potentials lower bounds     
% es_constraints.mu0_max            nm x 1  std chemical potentials upper bounds     
% es_constraints.mu_fix             nm x 1  chemical potentials predefined values
% es_constraints.mu_min             nm x 1  chemical potentials lower bounds     
% es_constraints.mu_max             nm x 1  chemical potentials upper bounds     
% es_constraints.dmu_fix            nr x 1  Reaction Gibbs free energies predefined values
% es_constraints.dmu_min            nr x 1  Reaction Gibbs free energies lower bounds  
% es_constraints.dmu_max            nr x 1  Reaction Gibbs free energies upper bounds  
% es_constraints.dmu_limit           1 x 1  Reaction Gibbs free energies upper limit for absolute values
% es_constraints.dmu_limit_min       1 x 1  Reaction Gibbs free energies lower limit for absolute values
% es_constraints.log_u_mean         nr x 1  ln enzyme levels data
% es_constraints.log_u_std          nr x 1  ln enzyme levels data error bars
% es_constraints.beta_M_fix         nr x nm saturation values beta_M predefined values
% es_constraints.beta_A_fix         nr x nm saturation values beta_A predefined values
% es_constraints.beta_I_fix         nr x nm saturation values beta_I predefined values
% es_constraints.alpha_min          (default 0.05) minimum value for uniform sampling of alpha values
% es_constraints.alpha_max          (default 0.95) maximum value for uniform sampling of alpha values
% es_constraints.alpha_A_mean       nr x nm mean values for alpha_A    (only for beta distribution)
% es_constraints.alpha_A_std        nr x nm std dev values for alpha_A (only for beta distribution)
% es_constraints.alpha_I_mean       nr x nm mean values for alpha_I    (only for beta distribution)
% es_constraints.alpha_I_std        nr x nm std dev values for alpha_I (only for beta distribution)
% es_constraints.alpha_M_mean       nr x nm mean values for alpha_M    (only for beta distribution)
% es_constraints.alpha_M_std        nr x nm std dev values for alpha_M (only for beta distribution)
% es_constraints.rho                 1 x 1  upper bound on v+/v or v-/v 
% es_constraints.dmu_eqconstraint  struct  Fields .matrix, .vector define equality constraints
% es_constraints.mu_eqconstraint   struct  Fields .matrix, .vector define equality constraints
%
% with nr = #reactions; nm = #metabolites; nx = #external metabolites
% For defining only some entries in a vector, replace the others by 'nan'
%
%
% 2. Structure 'es_options': general directives for the algorithm
%
% Name                        Type/default   Meaning 
% ----                        ------------   -------
% es_options.seed                            Random seed
% es_options.sampling_method                 Strategy for steady state sampling 
%                                              'v from data'
%                                              'v and mu'
%                                              'c0 and c'
% es_options.method_flux_sampling            Strategy for flux sampling 
%                                              'sample_and_discard'
%                                              'convex_optimisation'
%                                              'cycle_correction'
% es_options.n_flux_samples       1          Number of flux samples
% es_options.flux_units                      'dimensionless', 'mM/s', 'mol/s'
% es_options.n_samples            1          Number of samples in multiple sampling (model ensemble)
% es_options.n_dmu_samples        1          For each of these: number of reaction affinity samples
% es_options.n_saturation_samples 1          For each of these: number of elasticity samples
% es_options.kinetic_law          string     Type of modular rate law (see modular rate laws, e.g., 'cs', 'ms')
% es_options.h                    nr x 1     Predefined reaction cooperativities
% es_options.set_alpha_to_half    boolean    Set all saturation values to 0.5 instead of random sampling?
% es_options.set_alpha_nonuniform boolean    Draw alpha values from beta distribution?
%                                              Use distribution parameters from the sparse matrices
%                                              es_constraints.alpha_A_mean etc  (only nonzero entries are used):
% es_options.KV_prior_mean       = 1;        kV value to be substituted in inactive reactions
% es_options.limit_cooperativity = 2;        threshold value for reaction cooperativity
% es_options.no_equilibrium       boolean    If rates vanish, assume that zero enzyme levels are the reason
% es_options.ind_ignore           []         Indices of reactions to be ignored when computing thermodynamic 
%                                              loops (only in for flux correction by loop substraction)
% es_options.zc                   []         Metabolite derivative of objective function
% es_options.zv                   []         Flux derivative of objective function
% es_options_in_silico.limit_elasticities_in_multireactions = 1; normalise scaled elasticities in multireactions (eg biomass reaction) to an sum of absolute values = 1
% es_options.verbose              (boolean)  Verbose output?
% es_options.graphics_flag        (boolean)  Show graphics?
% es_options.print_graphics       (boolean)  Save graphics to file?
% es_options.flag_test            (boolean)  Run tests? (default 0)

% TO BE IMPLEMENTED (matrices should rather be in 'constraints')
% es_options.sample_K     Sample KM, KA, KI values instead if saturation values
%                         Use distribution parameters from the sparse matrices:
%                         (only nonzero entries are used):
% es_constraints.KA_mean  mean values for alpha_A
% es_constraints.KA_std   std dev values for alpha_A
% es_constraints.KI_mean  mean values for alpha_I
% es_constraints.KI_std   std dev values for alpha_I
% es_constraints.KM_mean  mean values for alpha_M
% es_constraints.KM_std   std dev values for alpha_M

if isstruct(N), N = N.N; end

[nm,nr] = size(N);

es_constraints.v_min        = - inf(nr,1);
es_constraints.v_max        =   inf(nr,1);
es_constraints.v_sign       = nan * ones(nr,1);
es_constraints.v_fix        = nan * ones(nr,1);
es_constraints.v_mean       = nan * ones(nr,1);
es_constraints.v_std        = nan * ones(nr,1);
es_constraints.ext_signs    = nan * ones(nm,1);

es_constraints.log_u_mean  = zeros(nr,1);
es_constraints.log_u_std   = zeros(nr,1);
es_constraints.log_c_mean  = zeros(nm,1);
es_constraints.log_c_std   = zeros(nm,1);
es_constraints.log_c_fix   = nan * ones(nm,1);
es_constraints.log_c_min   = nan * ones(nm,1);
es_constraints.log_c_max   = nan * ones(nm,1);
es_constraints.log_c0_mean = zeros(nm,1);
es_constraints.log_c0_std  = zeros(nm,1);

es_constraints.mu0_fix     = nan * ones(nm,1);
es_constraints.mu0_min     = nan * ones(nm,1);
es_constraints.mu0_max     = nan * ones(nm,1);

es_constraints.dmu_limit   = 100;
es_constraints.dmu_limit_min = 0.099;
es_constraints.dmu_min     = -es_constraints.dmu_limit*ones(nr,1);
es_constraints.dmu_max     =  es_constraints.dmu_limit*ones(nr,1);
es_constraints.dmu_fix     = nan*ones(nr,1);
es_constraints.dmu_eqconstraint = [];
es_constraints.mu_min      = -inf*ones(nm,1);
es_constraints.mu_max      =  inf*ones(nm,1);
es_constraints.mu_fix      = nan*ones(nm,1);
es_constraints.mu_eqconstraint = [];
es_constraints.Keq_fix     = nan*ones(nr,1);
es_constraints.beta_M_fix  = nan*ones(nr,nm);
es_constraints.beta_A_fix  = nan*ones(nr,nm);
es_constraints.beta_I_fix  = nan*ones(nr,nm);

es_constraints.alpha_min            = 0.05;
es_constraints.alpha_max            = 0.95;
es_constraints.alpha_A_mean         = [];
es_constraints.alpha_A_std          = [];   
es_constraints.alpha_I_mean         = [];
es_constraints.alpha_I_std          = [];
es_constraints.alpha_M_mean         = []; 
es_constraints.alpha_M_std          = [];

es_constraints.rho             = 100;
                               
es_options.seed                = nan;
es_options.sampling_method     = 'v and mu';
es_options.method_flux_sampling= 'sample_and_discard';
es_options.sampling_cycle_correction = 0;
es_options.n_flux_samples      = 1;
es_options.n_samples           = 10;
es_options.flux_units          = 'dimensionless';
es_options.n_dmu_samples       = 1;
es_options.n_saturation_samples= 1;
es_options.set_alpha_to_half   = 1;
es_options.set_alpha_nonuniform = 0;
es_options.kinetic_law          = 'cs';
es_options.verbose              = 0;
es_options.epsilon_stationary   = 10^-4;
es_options.cycles               = nan;
es_options.KV_prior_mean        = 1;  % value to be used where enzyme is switched off
es_options.h                    = ones(nr,1); % reaction cooperativity
es_options.limit_cooperativity  = 2;
es_options.zc                   = []; % metabolite derivative of target function
es_options.zv                   = []; % flux derivative of target function
es_options_in_silico.limit_elasticities_in_multireactions = 1; 
es_options.flag_second_order    = 1;
es_options.no_equilibrium       = 1;  
es_options.ind_ignore           = [];
es_options.graphics_flag        = 1;   
es_options.print_graphics       = 0;
es_options.flag_test            = 0;   
