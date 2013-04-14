function [r,r_orig,kinetic_data] = parameter_balancing_kinetic(network, kinetic_data, filename, Keq_given);

% [r,r_orig,kinetic_data] = parameter_balancing_kinetic(network, kinetic_data, filename, Keq_given);
%
% Convenience function for using parameter balancing to obtain 
% sensible kinetic constants, possibly given the equilibrium constants of a system
% 
% The code assumes that the network structure contains KEGG IDs and uses
% SBtab data files supplied that can contain data on 'standard chemical potential',
% 'Michaelis constant', 'activation constant',  'inhibitory constant','equilibrium constant',
% 'substrate catalytic rate constant', 'product catalytic rate constant'
% or some of these
%
% The standard reaction directions in the model have to follow the convention in KEGG!!!!

eval(default('filename','[]', 'Keq_given','[]'));

if 0,

  %% ---------
  %% example

  model_name = 'ecoli_glycolysis';
  filenames  = pathway_modelling_filenames(model_name);
  load(filenames.network_file);
  filename = {'~/projekte/pathway_modelling/data/GFEformation_flamholz_2011_millimolar.tsv', ...
              '~/projekte/pathway_modelling/data/kinetic_constants_brenda.tsv'};
  kinetic_data = [];
  
  %% kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','Michaelis constant', 'activation constant',  'inhibitory constant','equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'}, [], network,  filename, 0, 1);
  %filename = [];
  Keq_given = [];

 [r,r_orig,kinetic_data] = parameter_balancing_kinetic(network, kinetic_data, filename, Keq_given);

  gp               = struct('arrowsize',0.1,'actstyle','none','actprintnames',0,'flag_edges',1);
  gp.actstyle      ='none'; 
  gp.arrowsize     = 0.05;
  gp.arrowstyle='none'; gp.actstyle = 'fixed';

  figure(1); clf; netgraph_concentrations(network_CoHid,[],r.Keq,1,gp); title('Balanced: Equilibrium constants');
  figure(2); clf; netgraph_concentrations(network_CoHid,r.mu0,[],1,gp); title('Balanced: Standard chemical potentials');
  figure(3); clf; netgraph_concentrations(network_CoHid,[],r.KV,1,gp); title('Balanced: Velocity constants');
  gp.edgevalues = r.KM';  gp.edgestyle='fixed'; 
  figure(4); clf; netgraph_concentrations(network_CoHid,[],[],1,gp); title('Balanced: Michaelis constants');

  gp.edgevalues = [];; 
  figure(11); clf; netgraph_concentrations(network_CoHid,[],r_orig.Keq,1,gp); title('Data: Equilibrium constants');
  figure(12); clf; netgraph_concentrations(network_CoHid,r_orig.mu0,[],1,gp); title('Data: Standard chemical potentials');
  gp.edgevalues = r_orig.KM';  gp.edgestyle='fixed'; 
  figure(14); clf; netgraph_concentrations(network_CoHid,[],[],1,gp); title('Data: Michaelis constants');

end

if isempty(kinetic_data),
  kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','Michaelis constant', ...
                      'activation constant',  'inhibitory constant','equilibrium constant',...
                      'substrate catalytic rate constant', 'product catalytic rate constant'}, ...
                                                    [], network,  filename, 0, 1);
end  

[nm,nr] = size(network.N);


% ------------------------------------------------------------------------
% Determine consistent parameter set by parameter balancing

quantity_info     = data_integration_load_quantity_info;
model_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'
}';
data_quantities   = {'standard chemical potential','Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'
}';
basic_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant'}';

my_kinetic_data = data_integration_bounds_pseudovalues(kinetic_data,quantity_info,1,network);

if length(Keq_given),
  ind_finite = find(isfinite(Keq_given));
  my_kinetic_data.Keq.median(ind_finite) = Keq_given(ind_finite);
  my_kinetic_data.Keq.mean(ind_finite)   = Keq_given(ind_finite);
  my_kinetic_data.Keq.std(ind_finite)    = 0.001 * Keq_given(ind_finite);
  my_kinetic_data.Keq.mean_ln(ind_finite) = log(Keq_given(ind_finite));
  my_kinetic_data.Keq.std_ln  = log(1+  my_kinetic_data.Keq.std(ind_finite));
end

network.kinetics = set_kinetics(network, 'cs');
task = parameter_balancing_task(network, my_kinetic_data, quantity_info, model_quantities, basic_quantities);
res  = parameter_balancing(task, quantity_info, struct('insert_pseudo_values',0));

r.mu0   = res.kinetics_posterior_mode.mu0  ;
r.KV    = res.kinetics_posterior_mode.KV   ;
r.KM    = res.kinetics_posterior_mode.KM   ;
r.KA    = res.kinetics_posterior_mode.KA   ;
r.KI    = res.kinetics_posterior_mode.KI   ;
r.Keq   = res.kinetics_posterior_mode.Keq  ;
r.Kcatf = res.kinetics_posterior_mode.Kcatf;
r.Kcatr = res.kinetics_posterior_mode.Kcatr;

r_orig.mu0   = kinetic_data.mu0.median  ;
r_orig.KM    = kinetic_data.KM.median   ;
r_orig.KA    = kinetic_data.KA.median   ;
r_orig.KI    = kinetic_data.KI.median   ;
r_orig.Keq   = kinetic_data.Keq.median  ;
if length(Keq_given),
  r_orig.Keq(ind_finite) = Keq_given(ind_finite);
end
r_orig.Kcatf = kinetic_data.Kcatf.median;
r_orig.Kcatr = kinetic_data.Kcatr.median;