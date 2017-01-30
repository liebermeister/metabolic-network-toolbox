function [c, mu0, Keq, A, my_kinetic_data] = parameter_balancing_thermodynamic(network, v, kinetic_data_file, options)

% [c, v, mu0, Keq, A, my_kinetic_data] = parameter_balancing_thermodynamic(network, v, kinetic_data_file, options)
%
% Compute a thermodynamically feasible metabolic state for a given flux distribution
%
% This is a wrapper function that calls parameter_balancing_task and parameter_balancing
%
% For an alternative function (with all numerical data given as explicit function arguments;
% no data structure file), see 'thermo_pb'
% 
% The flux distribution must be thermodynamically feasible
% 
% The code assumes that the network structure contains metabolite KEGG IDs.
% It uses SBtab data files that can contain data on standard chemical potentials,
% equilibrium constants, concentrations, and reaction affinities
%
% The argument 'kinetic_data_file' can contain
%   o the name of a file with nuemrical input data, 
%   o a list of such kinetic_data_files
%   o OR: a kinetic data object, previously obtained from such a file using the function 'data_integration_load_kinetic_data.m'
%
% Fields of 'options':
%   options.ind_water           (indices of metabolites representing water)
%   options.set_water_conc_to_one = 1;
%   options.data_refer_to_molar = 0;    (flag)
%   options.flag_pseudo_value   = 0;    (flag; use pseudo values?)
%   options.A_max               = 1000; (kJ/mol)
%   options.A_min               = 0.5;  (kJ/mol)
%   options.A_mean              = nan;  (kJ/mol)
%   options.A_std               = std;  (kJ/mol)
%   options.A_lower             = nan * ones(nr,1);
%   options.A_upper             = nan * ones(nr,1);
%   options.conc_min            = 0.00001; (mM)
%   options.conc_max            = 100;     (mM)
%   options.c_min               = nan   (vector)
%   options.c_max               = nan   (vector)
%   options.c_fix               = nan   (vector)
%   options.c_fix_strict        = 0;    strictly fix concentrations (otherwise, allow for a range)!
%   options.u_max               = 1000; (mM)
%   options.u_min               = 0.01; (mM)
%   options.variability         = 2;    variability of known concentrations
%   options.sigma_mu0           = 3;    error of mu0 values (kJ/mol); 3 for alberty data
%   options_default.quantity_info_filename = []; file containing the prior table to be used
%   options_default.test_v_for_being_feasible = 1; run previous test for feasible v
%   options_default.fix_Keq_in_sampling = 0;
%
% my_kinetic_data: data and assumptions finally used in the optimisation task

[nm,nr] = size(network.N);

eval(default('options','struct'));

options_default.ind_water             = [];
options_default.set_water_conc_to_one = 1;
options_default.data_refer_to_molar   = 0;
options_default.flag_pseudo_values    = 0;
options_default.sigma_mu0             = 3;%    error of mu0 values (kJ/mol); 3 for alberty data
options_default.A_max                 = 100;
options_default.A_min                 = 0.5;  
options_default.c_min                 = nan * ones(nm,1);
options_default.c_max                 = nan * ones(nm,1);
options_default.c_fix                 = nan * ones(nm,1);
options_default.c_fix_strict          = 0;
options_default.A_fix                 = nan * ones(nr,1);
options_default.A_lower               = nan * ones(nr,1);
options_default.A_upper               = nan * ones(nr,1);
options_default.conc_min              = 0.00001; %(mM)
options_default.conc_max              = 100;     %(mM)
options_default.virtual_reactions     = {};
options_default.quantity_info_filename = [];
options_default.test_v_for_being_feasible = 1;
options_default.fix_Keq_in_sampling   = 0;

options = join_struct(options_default,options);

% --------------------------------------------------------

if options.test_v_for_being_feasible,
  %% dmu_abs_min    = 0;
  %% dmu_abs_max    = 10^10;
  %% [feasible_v, dmu, dmu_abs_max] = eba_feasible_lp(v,network.N,[],dmu_abs_min,[],dmu_abs_max);
  
  dmu_abs_min    = options.A_min;
  dmu_abs_max    = options.A_max;
  [feasible_v, dmu, dmu_abs_max] = eba_feasible_lp(v,network.N,[],dmu_abs_min,[],dmu_abs_max);
  if feasible_v==0,
    error('Infeasible flux distribution at given thermodynamic constraints');
  end
  
end

% --------------------------------------------------------
% virtual reactions, whose affinities should be controlled
% (e.g., the driving force in spontaneous ATP hydrolysis)

if length(options.virtual_reactions),

  if isstruct(kinetic_data_file),
    error('If virtual reactions are used, the kinetic data cannot be predefined, but must be given in the form of files.');
  end
    
  %% add virtual reactions to the network and run parameter balancing 
  v_aug           = [v; ones(length(options.virtual_reactions),1)];
  options.A_lower = [options.A_lower; options.virtual_A_lower];
  options.A_upper = [options.A_upper; options.virtual_A_upper];

  network_aug     = network;
  for it = 1:length(options.virtual_reactions),
    nn = zeros(nm,1);
    ind     = label_names(options.virtual_reactions(it).metabolites, network.metabolites);
    nn(ind) = options.virtual_reactions(it).stoichiometries;
    network_aug.N       = [network_aug.N, nn];
    network_aug.actions = [network_aug.actions; {sprintf('aug_%f',it)}];
  end

  ooptions = rmfield(options, 'virtual_reactions');
  [c, mu0, Keq_aug, A_aug, my_kinetic_data] = parameter_balancing_thermodynamic(network_aug, v_aug, kinetic_data_file, ooptions);

  A   =   A_aug(1:end-length(options.virtual_reactions));
  Keq = Keq_aug(1:end-length(options.virtual_reactions));

  return

end


% ------------------------------------------------------------------------
% read standard chemical potentials 

% load Gibbs free energies of formation ... this requires metabolite KEGG IDs in the model!

% the variable 'kinetic_data_file' itelf can already contain the data: 

if isstruct(kinetic_data_file),
  if isfield(kinetic_data_file,'mu0'), kinetic_data.mu0 = kinetic_data_file.mu0; end 
  if isfield(kinetic_data_file,'Keq'), kinetic_data.Keq = kinetic_data_file.Keq; end 
  if isfield(kinetic_data_file,'c'),   kinetic_data.c = kinetic_data_file.c; end 
  if isfield(kinetic_data_file,'A'),   kinetic_data.A = kinetic_data_file.A; end 
else
  data_quantities = {'standard chemical potential', 'equilibrium constant', 'concentration', 'reaction affinity'}';
  quantity_info   = data_integration_load_quantity_info;
  kinetic_data    = data_integration_load_kinetic_data(data_quantities, quantity_info, network, kinetic_data_file, 0, 1);
end

% for which metabolite concentrations do we have standard chemical potentials?
% mu0 = kinetic_data.mu0.median;
% figure(2); netgraph_concentrations(network_CoHid,isfinite(mu0),[],1,gp) 

if options.data_refer_to_molar,

  ind_cor = setdiff(1:nm,options.ind_water);
  kinetic_data.mu0.median(ind_cor)  = kinetic_data.mu0.median(ind_cor) - RT * log(1000); 
  kinetic_data.mu0.mean(ind_cor)    = kinetic_data.mu0.mean(ind_cor)   - RT * log(1000); 
  kinetic_data.mu0.lower(ind_cor)   = kinetic_data.mu0.lower(ind_cor)  - RT * log(1000); 
  kinetic_data.mu0.upper(ind_cor)   = kinetic_data.mu0.upper(ind_cor)  - RT * log(1000); 
  kinetic_data.c.median             = kinetic_data.c.median   * 1000; 
  kinetic_data.c.mean               = kinetic_data.c.mean     * 1000; 
  kinetic_data.c.std                = kinetic_data.c.std      * 1000; 
  kinetic_data.c.lower              = kinetic_data.c.lower    * 1000; 
  kinetic_data.c.upper              = kinetic_data.c.upper    * 1000; 
  kinetic_data.c.mean_ln            = kinetic_data.c.mean_ln  + log(1000); 
  kinetic_data.c.lower_ln           = kinetic_data.c.lower_ln + log(1000); 
  kinetic_data.c.upper_ln           = kinetic_data.c.upper_ln + log(1000); 
  
  kinetic_data = concentration_sampling_update_kinetic_data(kinetic_data, network, v, options);

end

if options.set_water_conc_to_one, 
  options.c_fix(options.ind_water) = 1; 
end

% -----------------------
% fix mu0 values and concentrations exactly at data values

% TOO SMALL VALUES HERE CAN LEAD TO NUMERICAL PROBLEMS LATER ON 
% IN PARAMETER BALANCING (QUADPROG)

kinetic_data.mu0.std(find(isfinite(kinetic_data.mu0.std)))   = 10^-4;

if isfield(kinetic_data,'c'),
  kinetic_data.c.std(find(isfinite(kinetic_data.c.std)))       = 10^-4;
  kinetic_data.c.std_ln(find(isfinite(kinetic_data.c.std_ln))) = 10^-4;
end


% ------------------------------------------------------------------------
% Determine consistent parameter set by parameter balancing

quantity_info     = data_integration_load_quantity_info([],options.quantity_info_filename);

model_quantities  = {'standard chemical potential','equilibrium constant', 'concentration','reaction affinity'}';
data_quantities   = {'standard chemical potential','equilibrium constant', 'concentration','reaction affinity'}';
basic_quantities  = {'standard chemical potential','concentration'}';

my_kinetic_data = kinetic_data;


% prepare data for parameter balancing

quantity_info.LowerBound{quantity_info.symbol_index.c}   = '0.001';
quantity_info.UpperBound{quantity_info.symbol_index.c}   = '100';
quantity_info.LowerBound{quantity_info.symbol_index.A}   = '-60';
quantity_info.UpperBound{quantity_info.symbol_index.A}   = '60';
quantity_info.LowerBound{quantity_info.symbol_index.Keq} = '0.000001';
quantity_info.UpperBound{quantity_info.symbol_index.Keq} = '1000000';
quantity_info.PriorStd{quantity_info.symbol_index.Keq}   = '0.05';

my_kinetic_data.mu0.std = options.sigma_mu0 * ones(nm,1);

my_kinetic_data = data_integration_bounds_pseudovalues(my_kinetic_data,quantity_info,options.flag_pseudo_values,network);

options.c_min(isnan(options.c_min)) = options.conc_min;
options.c_max(isnan(options.c_max)) = options.conc_max;

my_kinetic_data.c.lower = options.c_min;
my_kinetic_data.c.upper = options.c_max;

ind_fix = find(isfinite(options.c_fix));
if length(ind_fix),
  if options.c_fix_strict, 
    my_kinetic_data.c.lower(ind_fix) = options.c_fix(ind_fix);
    my_kinetic_data.c.upper(ind_fix) = options.c_fix(ind_fix);
  else
    display('Allowing for variation around fixed concentrations in parameter_balancing_thermodynamic.m - to change this, set c_fix_strict = 1;');
    my_kinetic_data.c.lower(ind_fix) = 0.95 * options.c_fix(ind_fix);
    my_kinetic_data.c.upper(ind_fix) = 1.05 * options.c_fix(ind_fix);
  end
end

my_kinetic_data.c.lower_ln = log(my_kinetic_data.c.lower);
my_kinetic_data.c.upper_ln = log(my_kinetic_data.c.upper);

my_kinetic_data.A.upper(v<0) = min(my_kinetic_data.A.upper(v<0),  -options.A_min);
my_kinetic_data.A.lower(v>0) = max(my_kinetic_data.A.lower(v>0),   options.A_min);
my_kinetic_data.A.upper(v>0) = max(my_kinetic_data.A.median(v>0),  options.A_max);
my_kinetic_data.A.lower(v<0) = min(my_kinetic_data.A.median(v<0), -options.A_max);

ind = find(isfinite(options.A_lower));
my_kinetic_data.A.lower(ind) = max(my_kinetic_data.A.lower(ind), options.A_lower(ind));

ind = find(isfinite(options.A_upper));
my_kinetic_data.A.upper(ind) = min(my_kinetic_data.A.upper(ind), options.A_upper(ind));

ind_plus = find(options.A_fix>0);
my_kinetic_data.A.lower(ind_plus) = 0.99 * options.A_fix(ind_plus);
my_kinetic_data.A.upper(ind_plus) = 1.01 * options.A_fix(ind_plus);

ind_minus = find(options.A_fix<0);
my_kinetic_data.A.lower(ind_minus) = 1.01 * options.A_fix(ind_minus);
my_kinetic_data.A.upper(ind_minus) = 0.99 * options.A_fix(ind_minus);

network.kinetics = set_kinetics(network, 'cs');


% -----------------------------------------------------
% Parameter balancing calculation

% display the adjusted data
% data_integration_display_kinetic_data(my_kinetic_data,network);

task   = parameter_balancing_task(network, my_kinetic_data, quantity_info, model_quantities, basic_quantities);

result = parameter_balancing(task, quantity_info, struct('use_pseudo_values',0,'fix_Keq_in_sampling',options.fix_Keq_in_sampling));


% -----------------------------------------------------
% Output variables

mu0 = result.kinetics_posterior_mode.mu0;
Keq = result.kinetics_posterior_mode.Keq;
c   = result.kinetics_posterior_mode.c;    
A   = result.kinetics_posterior_mode.A;


% -----------------------------------------------------
% Test whether results are consistent

if norm( [sign(network.N' *[mu0 + RT * log(c)]) + sign(v)] .* double([v~=0])),
  task.network.actions(find([sign(network.N' *[mu0 + RT * log(c)]) + sign(v)]))
  inconsistent_solution = 1;
else
  inconsistent_solution = 0;
end

if norm(network.N'*mu0/RT + log(Keq)) / norm(log(Keq)) > 10^-5,
  warning('WEGSCHEIDER CONDITION IS NOT EXACTLY SATISFIED')
  if inconsistent_solution == 0,
    display('But solution is consistent with flux directions');
  end
end

if norm([network.N' *[mu0 + RT * log(c)] + A]) / norm(A) > 10^-5,
  warning('RELATIONSHIP BETWEEN C, KEQ AND A IS NOT EXACTLY SATISFIED')
  %[network.N' *[mu0 + RT * log(c)], - A]
  if inconsistent_solution == 0,
    display('But solution is consistent with flux directions');
  end
end

if inconsistent_solution,
  error('Solution found is inconsistent with flux directions');
end