function [options, constraints] = flux_constraints_default_options(N)

% FLUX_CONSTRAINTS_DEFAULT_OPTIONS - Default flux and thermodynamic constraints and general options
%
% [options, constraints] = flux_constraints_default_options(N)
%
% Initialise the structures 'options' and 'constraints' used in flux sampling and thermodynamic
% calculations (e.g. sample_feasible_v, feasible_fluxes, sample_feasible_dmu, mnt_mu_simple_guess,
% mnt_check_flux_constraints). Structural thermokinetic modelling (STM, function stm_default_options
% in the STM toolbox) extends these structures by further fields.
%
% The default values are read from the SBtab Config file 'flux-constraints_Config.tsv' (in the same folder
% as this function; TableID 'ConfigureFluxConstraints'; option IDs 'options.<name>' and 'constraints.<name>').
% To change the defaults, edit this file. For vectors, the Config file contains a single value that is used
% for all elements; the default reaction Gibbs free energy bounds dmu_min and dmu_max (empty in the Config
% file) are given by -dmu_limit and dmu_limit.
%
% Function arguments:
% N can be the (complete) stoichiometric matrix or a 'network' data structure
%
% Function outputs:
%
% 1. Structure 'constraints'
%
% Field name                    Size    Meaning
% ----------                    ----    -------
% constraints.v_fix              nr x 1  fluxes predefined values (overrides vmin and vmax)
% constraints.v_min              nr x 1  fluxes lower bounds (overrides -vmax)
% constraints.v_max              nr x 1  fluxes upper bounds
% constraints.v_sign             nr x 1  fluxes signs (overrides vmin and vmax)
% constraints.v_mean             nr x 1  fluxes data (can also be a matrix)
% constraints.v_std              nr x 1  fluxes data error bars (can also be a matrix)
% constraints.Keq_fix            nr x 1  equilibrium constants predefined values
% constraints.ext_signs          nm x 1  signs of production rates of external metabolites (nan for internal metabolites)
% constraints.log_c_fix          nm x 1  ln metabolite levels predefined values
% constraints.log_c_min          nm x 1  ln metabolite levels lower bounds
% constraints.log_c_max          nm x 1  ln metabolite levels upper bounds
% constraints.log_c_mean         nm x 1  ln metabolite levels data
% constraints.log_c_std          nm x 1  ln metabolite levels data error bars
% constraints.mu0_fix            nm x 1  std chemical potentials predefined values
% constraints.mu0_min            nm x 1  std chemical potentials lower bounds
% constraints.mu0_max            nm x 1  std chemical potentials upper bounds
% constraints.mu_fix             nm x 1  chemical potentials predefined values
% constraints.mu_min             nm x 1  chemical potentials lower bounds
% constraints.mu_max             nm x 1  chemical potentials upper bounds
% constraints.dmu_fix            nr x 1  Reaction Gibbs free energies predefined values
% constraints.dmu_min            nr x 1  Reaction Gibbs free energies lower bounds
% constraints.dmu_max            nr x 1  Reaction Gibbs free energies upper bounds
% constraints.dmu_limit           1 x 1  Reaction Gibbs free energies upper limit for absolute values
%                                        (when raising it, raise dmu_min and dmu_max as well)
% constraints.dmu_limit_min       1 x 1  Reaction Gibbs free energies lower limit for absolute values
% constraints.rho                 1 x 1  upper bound on v+/v or v-/v
% constraints.dmu_eqconstraint  struct  Fields .matrix, .vector define equality constraints
% constraints.mu_eqconstraint   struct  Fields .matrix, .vector define equality constraints
%
% with nr = #reactions; nm = #metabolites
% For defining only some entries in a vector, replace the others by 'nan'
%
% 2. Structure 'options'
%
% Name                        Type      Meaning
% ----                        ----      -------
% options.seed                          Random seed (nan: no seed)
% options.verbose             boolean   Verbose output?
% options.epsilon_stationary  float     Tolerance for stationarity
% options.cycles                        Thermodynamic cycles (NaN: computed when needed)
% options.n_dmu_samples       integer   Number of reaction Gibbs free energy samples (in sample_feasible_dmu)
% options.ind_ignore          vector    Indices of reactions to be ignored when computing thermodynamic
%                                        loops (only in for flux correction by loop substraction)

if isstruct(N), N = N.N; end

[nm,nr] = size(N);

config_file = fullfile(fileparts(mfilename('fullpath')), 'flux-constraints_Config.tsv');
config      = load_config_table(config_file);
options     = config.options;
constraints = config.constraints;

% expand the values from the Config file to vectors

for f = {'v_min','v_max','v_sign','v_fix','v_mean','v_std','dmu_fix','Keq_fix'}
  constraints.(f{1}) = constraints.(f{1}) * ones(nr,1);
end

for f = {'ext_signs','log_c_mean','log_c_std','log_c_fix','log_c_min','log_c_max', ...
         'mu0_fix','mu0_min','mu0_max','mu_min','mu_max','mu_fix'}
  constraints.(f{1}) = constraints.(f{1}) * ones(nm,1);
end

if isempty(constraints.dmu_min), constraints.dmu_min = -constraints.dmu_limit; end
if isempty(constraints.dmu_max), constraints.dmu_max =  constraints.dmu_limit; end
constraints.dmu_min = constraints.dmu_min * ones(nr,1);
constraints.dmu_max = constraints.dmu_max * ones(nr,1);
