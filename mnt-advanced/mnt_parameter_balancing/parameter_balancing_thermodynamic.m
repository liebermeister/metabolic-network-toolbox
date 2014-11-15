function [c, mu0, Keq, A, my_kinetic_data] = parameter_balancing_thermodynamic(network, v, filename, options)

% [c, v, mu0, Keq, A, my_kinetic_data] = parameter_balancing_theromdynamic(network, v, filename,options)
%
% Convenience function for using parameter balancing to obtain a steady state 
% given a thermodynamically feasible flux distribution
% 
% The code assumes that the network structure contains KEGG IDs and uses
% SBtab data files supplied that can contain data on standard chemical potentials,
% equilibrium constants, concentrations, and reaction affinities
%
% Fields of options:
%   options.ind_water           (indices of metabolites representing water)
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
%   options.c_fix               = nan   (vector)
%   options.u_max               = 1000; (mM)
%   options.u_min               = 0.01; (mM)
%   options.variability         = 2;    variability of known concentrations
%   options.sigma_mu0           = 3;    error of mu0 values (kJ/mol); 3 for alberty data
%
% my_kinetic_data: data and assumptions finally used in the optimisation task

[nm,nr] = size(network.N);

eval(default('options','struct'));

options_default.ind_water           = [];
options_default.set_water_conc_to_one = 1;
options_default.data_refer_to_molar = 0;
options_default.flag_pseudo_values  = 0;
options_default.sigma_mu0           = 3;%    error of mu0 values (kJ/mol); 3 for alberty data
options_default.A_max               = 100;
options_default.A_min               = 0.5;  
options_default.c_fix               = nan * ones(nm,1);
options_default.A_fix               = nan * ones(nr,1);
options_default.A_lower             = nan * ones(nr,1);
options_default.A_upper             = nan * ones(nr,1);
options_default.conc_min            = 0.00001; %(mM)
options_default.conc_max            = 100;     %(mM)
options_default.virtual_reactions   = {};

options = join_struct(options_default,options);


% --------------------------------------------------------
% virtual reactions, whose affinities should be controlled
% (e.g., the driving force in spontaneous ATP hydrolysis)

if length(options.virtual_reactions),

  %% add virtual reactions to the network and run parameter balancing 
  v_aug           = [v; ones(length(options.virtual_reactions),1)];
  options.A_lower = [options.A_lower; options.virtual_A_lower];
  options.A_upper = [options.A_upper; options.virtual_A_upper];

  network_aug     = network;
  for it = 1:length(options.virtual_reactions),
    nn = zeros(nm,1);
    ind     = label_names(options.virtual_reactions(it).metabolites, network.metabolites);
    nn(ind) = options.virtual_reactions(it).stoichiometries;
    network_aug.N   = [network_aug.N, nn];
    network_aug.actions = [network_aug.actions; {sprintf('aug_%f',it)}];
  end

  ooptions = rmfield(options, 'virtual_reactions');
  [c, mu0, Keq_aug, A_aug, my_kinetic_data] = parameter_balancing_thermodynamic(network_aug, v, filename, ooptions);

  A   =   A_aug(1:end-length(options.virtual_reactions));
  Keq = Keq_aug(1:end-length(options.virtual_reactions));

  return

end


% ------------------------------------------------------------------------
% read standard chemical potentials 

% load Gibbs free energies of formation ... this requires metabolite KEGG IDs in the model!

data_quantities = {'standard chemical potential','equilibrium constant', 'concentration','reaction affinity'}';
quantity_info   = data_integration_load_quantity_info;
kinetic_data    = data_integration_load_kinetic_data(data_quantities, quantity_info, network, filename, 0, 1);

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
kinetic_data.c.std(find(isfinite(kinetic_data.c.std)))       = 10^-4;
kinetic_data.c.std_ln(find(isfinite(kinetic_data.c.std_ln))) = 10^-4;


% ------------------------------------------------------------------------
% Determine consistent parameter set by parameter balancing

quantity_info     = data_integration_load_quantity_info;
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

my_kinetic_data.c.lower    = options.conc_min * ones(nm,1);
my_kinetic_data.c.upper    = options.conc_max * ones(nm,1);

ind_fix = find(isfinite(options.c_fix));
my_kinetic_data.c.lower(ind_fix) = 0.95 * options.c_fix(ind_fix);
my_kinetic_data.c.upper(ind_fix) = 1.05 * options.c_fix(ind_fix);
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

task = parameter_balancing_task(network, my_kinetic_data, quantity_info, model_quantities, basic_quantities);

res  = parameter_balancing(task, quantity_info, struct('insert_pseudo_values',0));


% -----------------------------------------------------
% Output variables

mu0 = res.kinetics_posterior_mode.mu0;
Keq = res.kinetics_posterior_mode.Keq;
c   = res.kinetics_posterior_mode.c;    
A   = res.kinetics_posterior_mode.A;
