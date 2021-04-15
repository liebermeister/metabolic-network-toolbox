function [network, r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std, r_geom_mean, r_geom_std, task, result] = parameter_balancing_sbtab(model_file, data_file, pb_options)

% PARAMETER_BALANCING_SBTAB - Wrapper function for parameter balancing 
%
%  [network, r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std, r_geom_mean, r_geom_std, task, result] = parameter_balancing_sbtab(model_file, data_file, pb_options)
%
% Wrapper function for parameter balancing, with filenames (for model and data) as input arguments
%  o reads input data (model, priors, kinetic and metabolic data) from SBtab files (the SBtab model file can be replaced by an SBML file)
%  o runs parameter balancing (without and with concentration data)
% 
% Function arguments:
%  model_file:           SBtab or SBML model filename; files with ".xml" extension are assumed to be SBML, otherwise SBtab
%  data_file:            SBtab data file (kinetic and other constants)
%  pb_options:           (optional) struct with settings for parameter balancing; for default values, see 'parameter_balancing_options'
%  pb_options.parameter_prior_file:        File with parameter priors (optional)
%  pb_options.use_python_version_defaults  Flag - use same default settings as in the python version
%
% Output data structures of data type 'kinetics' 
%  r            posterior mode (considering constraints)
%  r_orig       original values
%  r_mean       posterior mode (ignoring constraints, may be outside the feasible range)
%  r_std        posterior std  (per parameter; ignoring constraints, >= actual posterior std)
%  r_geom_mean  posterior mode (ignoring constraints, may be outside the feasible range)
%  r_geom_std   posterior geom std  (per param.; ignoring constraints, >= actual posterior geom std)
%  r_samples    cell array of data struct with sampled parameter values
%
% With default options, the function uses the information from input files; 
% Set the option 'pb_options.use_python_version_defaults=1' for a direct comparison to the python parameter balancing tool
%
% This function calls the functions 'parameter_balancing_task', 'parameter_balancing', and 'parameter_balancing_output'

eval(default('pb_options', 'struct'));

pb_options = parameter_balancing_update_options(join_struct(parameter_balancing_options, pb_options));

log_text = '';

% ----------------------------------------------------------
% load model

display(sprintf('o Using model file %s', model_file))

network = network_import_model(model_file);

if isfield(network,'metabolite_is_an_enzyme'),
  if sum(network.metabolite_is_an_enzyme),
    warning('At least one of the model species appears to be an enzyme. Parameter balancing currently cannot handle this case. Please remove all enzyme species from your model and run parameter balancing again.')
  end
end

% ----------------------------------------------------------
% load parameter_prior and define relevant quantities

parameter_prior = parameter_balancing_prior([],pb_options.parameter_prior_file,1); 
parameter_prior = pb_parameter_prior_adjust(parameter_prior, pb_options); 

[model_quantities, basic_quantities, data_quantities, pseudo_quantities] = parameter_balancing_quantities(parameter_prior, network, pb_options);


% ----------------------------------------------------------
% load kinetic data

if length(data_file), 
  if iscell(data_file),
    for it = 1:length(data_file),
      display(sprintf('o Using data file %s', data_file{it})); 
    end
  else
    display(sprintf('o Using data file %s', data_file)); 
  end
end

kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, data_file, struct('use_sbml_ids', pb_options.use_sbml_ids, 'use_kegg_ids', pb_options.use_kegg_ids, 'use_python_version_defaults', pb_options.use_python_version_defaults, 'verbose', pb_options.verbose));

kinetic_data_orig = kinetic_data;
kinetic_data      = pb_kinetic_data_adjust(kinetic_data, parameter_prior, network, pb_options);

% Display the adjusted data
% parameter_balancing_kinetic_data_show(kinetic_data);

% -----------------------------------------------------------
% Run parameter balancing

task = parameter_balancing_task(network, kinetic_data, parameter_prior, model_quantities, basic_quantities, pseudo_quantities);

result = parameter_balancing_calculation(task, parameter_prior, pb_options);

[r, r_mean, r_std, r_geom_mean, r_geom_std, r_orig, r_samples] = parameter_balancing_output(result, kinetic_data_orig, pb_options, network);

if pb_options.fix_irreversible_reactions,
  %% irreversible reactions are reactions whose log10 Keq value, divided by 
  %% the sum of substrate and product molecularities, is below -5 or above 5.
  %% for our test, we consider the data Keq values (if a value is nan, we replace it by the balanced Keq for the same reaction)
  molecularity_sums  = sum(abs(network.N),1)';
  log10_Keq          = log10(r_orig.Keq);
  log10_Keq(~isfinite(r_orig.Keq)) = log10(r_orig.Keq(~isfinite(r_orig.Keq)));
  ind_irreversible   = find(abs(log10_Keq) ./ molecularity_sums > 5);
  if length(ind_irreversible),
    sign_irreversible  = sign(log10_Keq(ind_irreversible));
    display('parameter_balancing_sbtab: setting backward Kcat values of supposedly irreversible reactions (and forward Kcat values of supposedly irreversible reactions in reverse direction) to 0');
    network.actions(ind_irreversible);
    r.Kcatr(sign_irreversible==1)  = 0;
    r.Kcatf(sign_irreversible==-1) = 0;
  end
end

network.kinetics = r;
