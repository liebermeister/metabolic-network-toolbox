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
%   options.A_min               = 0.01; (kJ/mol)
%   options.A_mean              = nan;  (kJ/mol)
%   options.A_std               = std;  (kJ/mol)
%   options.conc_min            = 0.;   (mM)
%   options.conc_max            = 10^10;(mM)
%   options.u_max               = 1000; (mM)
%   options.u_min               = 0.01; (mM)
%   options.variability         = 2;    variability of known concentrations
%   options.sigma_mu0           = 3;    error of mu0 values (kJ/mol); 3 for alberty data
%
% my_kinetic_data: data and assumptions finally used in the optimisation task
%
% Code originally from: hepatonet1 ...

eval(default('options','struct'));

options_default.data_refer_to_molar = 0;
options_default.flag_pseudo_values  = 0;
options_default.sigma_mu0           = 3;%    error of mu0 values (kJ/mol); 3 for alberty data
options_default.A_max               = 1000;
options_default.A_min               = 0.01;

options = join_struct(options_default,options);


% ------------------------------------------------------------------------

if 0,

  % ---------
  model_name     = 'eccm'; 
  condition_name = 'ishii_growth_07';
  filenames      = elasticity_sampling_filenames(model_name,condition_name); 
  cd([filenames.es_dir]);  
  load(filenames.model_data_file);
  load(filenames.reference_state_file);
  v = result.v;
  options.ind_water = find(strcmp('H2O',network.metabolites));
  filename        = {'/home/wolfram/projekte/pathway_modelling/data/GFEformation_flamholz_2011_millimolar.tsv', ...
                     '/home/wolfram/projekte/pathway_modelling/data/concentration_ecoli_glucose.tsv', ...
                     '/home/wolfram/projekte/pathway_modelling/models/ecoli_glycolysis2/ecoli_glycolysis2_QuantityData.tsv', ...
                     '/home/wolfram/projekte/pathway_modelling/data/ATPase_affinity.tsv'}; 
  options.data_refer_to_molar = 0;

  filename = '/home/wolfram/projekte/pathway_modelling/models/ecoli_glycolysis2/ecoli_glycolysis2_QuantityData.tsv';
  options.data_refer_to_molar = 1;

  options.conc_min    = 0.; % in mM
  options.conc_max    = 10^10;   % in mM
  options.A_max       = 1000;
  options.A_min       = 0.;   % minimal Delta G needed to drive a forward reaction
  
  % ---------
  model_name     = 'ecoli_glycolysis2'; 
  condition_name = 'glucose';
  filenames = pathway_modelling_filenames(model_name);
  load(filenames.network_file);
  load(filenames.flux_file);
  filename  = filenames.table_quantities;
  options.ind_water = find(strcmp('H2O',network.metabolites));
  options.data_refer_to_molar = 1;

  options.conc_min    = 0.001; % in mM
  options.conc_max    = 100;   % in mM
  options.A_max       = 60;
  options.A_min       = 0.5;   % minimal Delta G needed to drive a forward reaction

  % ---------
  options.variability = 2;     % variability of known concentrations
  options.sigma_mu0   = 3;     % error of mu0 values (kJ/mol); 3 for alberty data
  options.A_mean      = nan;   % mean and std values for forward directions
  options.A_std       = nan;

  [c, mu0, Keq, A, my_kinetic_data] = parameter_balancing_thermodynamic(network, v, filename,options);

  gp               = struct('arrowsize',0.1,'actstyle','none','actprintnames',0);
  gp.omitreactions = {'Biomass production'};
  gp.actstyle      ='none'; gp.arrowstyle='fluxes';
  gp.arrowsize     = 0.05;

  %figure(1); netgraph_concentrations(network_CoHid,[],sign(v),1,gp); title('Fluxes');
  figure(1); netgraph_concentrations(network_CoHid,[],v,1,gp); title('Fluxes');
  gp.arrowstyle='none'; gp.actstyle = 'fixed';
  figure(2); netgraph_concentrations(network_CoHid,log10(c),[],1,gp); title('log10 Concentrations');
  %figure(4); netgraph_concentrations(network_CoHid,mu0,[],1,gp); title('Standard chemical potentials'); 
  %figure(3); netgraph_concentrations(network_CoHid,mu0+RT*log(c),[],1,gp); title('Chemical potentials');
  gp.arrowvalues=A; gp.actstyle = 'fixed'; gp.arrowstyle = 'fluxes'; gp.arrowsize= 0.05;
  figure(5); netgraph_concentrations(network_CoHid,[],abs(A),1,gp); title('Reaction affinities'); 
  %% changes of concentrations during balancing
  figure(100); plot(my_kinetic_data.c.median,c,'.'); set(gca,'XScale','log','YScale','log');

end


% ------------------------------------------------------------------------

[nm,nr] = size(network.N);


% ------------------------------------------------------------------------
% read standard chemical potentials 

% load gibbs free energies of formation ... this requires metabolite KEGG IDs in the model!

data_quantities = {'standard chemical potential','equilibrium constant', 'concentration','reaction affinity'}';
quantity_info   = data_integration_load_quantity_info;
kinetic_data    = data_integration_load_kinetic_data(data_quantities, quantity_info, network, filename, 0, 1);

% for which metabolites do we have standard chemical potentials?
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

my_kinetic_data.c.lower    = 0.0001 * ones(nm,1);
my_kinetic_data.c.upper    = 100 * ones(nm,1);
my_kinetic_data.c.lower_ln = log(0.0001) * ones(nm,1);
my_kinetic_data.c.upper_ln = log(100) * ones(nm,1);

my_kinetic_data.A.upper(v<0) = min(my_kinetic_data.A.upper(v<0),-options.A_min);
my_kinetic_data.A.lower(v>0) = max(my_kinetic_data.A.lower(v>0), options.A_min);
my_kinetic_data.A.upper(v>0) = options.A_max;
my_kinetic_data.A.lower(v<0) = -options.A_max;

network.kinetics = set_kinetics(network, 'cs');
task = parameter_balancing_task(network, my_kinetic_data, quantity_info, model_quantities, basic_quantities);
res  = parameter_balancing(task, quantity_info, struct('insert_pseudo_values',0));

mu0 = res.kinetics_posterior_mode.mu0;
Keq = res.kinetics_posterior_mode.Keq;
c   = res.kinetics_posterior_mode.c;    
A   = res.kinetics_posterior_mode.A;
