function [network, r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std, r_geom_mean, r_geom_std] = parameter_balancing_sbtab(model_file, data_file, pb_options)

% PARAMETER_BALANCING_SBTAB - Wrapper function for parameter balancing 
%
%  [network, r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std] = parameter_balancing_sbtab(model_file, data_file, pb_options)
%
% With no options set, the function uses directly the information from input files; it is best suited for comparisons with other PB tools
%
%  o reads input data (model, priors, kinetic and metabolic data) from SBtab files
%  o runs parameter balancing (without and with concentration data)
% 
% Function arguments:
%  o model_file:                   SBtab or SBML model filename; files with ".xml" extension are assumed to be SBML, otherwise SBtab
%  o data_file:            SBtab data file (kinetic and other constants)
%  o pb_options.parameter_prior_file: File with parameter priors (optional)
% 
% pb_options: struct with options; for default values, see 'parameter_balancing_options'
%
% This function calls the functions 'parameter_balancing_task', 'parameter_balancing', and 'parameter_balancing_output'

pb_options = parameter_balancing_update_options(join_struct(parameter_balancing_options, pb_options));

% ----------------------------------------------------------
% load model
  
if strcmp(model_file(end-3:end), '.xml'),
  display(sprintf('o Using model file %s', model_file))
  network = network_sbml_import(model_file);
else,
  network = sbtab_to_network(model_file, struct('kinetic_law','cs'));
end

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

if length(data_file), display(sprintf('o Using data file %s', data_file)); end  

kinetic_data      = data_integration_load_kinetic_data(data_quantities, parameter_prior, network, data_file, pb_options.use_sbml_ids, pb_options.use_kegg_ids);
kinetic_data_orig = kinetic_data;
kinetic_data      = pb_kinetic_data_adjust(kinetic_data, parameter_prior, network, pb_options);
% Show kinetic data importedfile_kinetic_data, 
% parameter_balancing_kinetic_data_show(kinetic_data);


% -----------------------------------------------------------
% run simple parameter balancing without metabolic data

task = parameter_balancing_task(network, kinetic_data, parameter_prior, model_quantities, basic_quantities, pseudo_quantities);

res  = parameter_balancing_calculation(task, parameter_prior,pb_options);

[r,r_mean,r_std,r_geom_mean,r_geom_std,r_orig,r_samples] = parameter_balancing_output(res,kinetic_data_orig,pb_options);

network.kinetics            = r;
