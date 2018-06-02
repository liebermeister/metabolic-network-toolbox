function balanced_parameters_SBtab = parameter_balancing(model_file, output_file, data_file, prior_file, options_file, model_name, options)

% PARAMETER_BALANCING Read model from SBML file; read data and general options from SBtab files; write results to output files
%
% parameter_balancing(model_file, output_file, data_file, prior_file, options_file, model_name, options)
%
% Additional options can be provided as a matlab struct "options". 
% (for a list of possible options, see "parameter_balancing_options")
%
% Arguments:
%
% model_file:     Model in SBML format (path and filename)
% output_file:    File basename (including directory path) for output files
% data_file:      (optional) Data in SBtab format (path and filename)
% prior_file:     (optional) Prior table in SBtab format (path and filename)
% options_file:   (optional) Options table in SBtab format (path and filename)
% model_name:     (optional) Model name 
% options:        (optional) Matlab struct with options, overriding options given in [options_file]
%  
% Optional arguments can be left empty ('[]')

% ----------------------------------------------------------------------
% Options
  
eval(default('data_file', '[]', 'prior_file', '[]', 'options_file', '[]', 'model_name', '[]', 'options','struct'));

options_default   = join_struct(parameter_balancing_options, struct('flag_check', 0));
options_from_file = sbtab_table_convert_to_simple_struct(sbtab_table_load(options_file),'Option','Value'); 
options           = join_struct(join_struct(options_default,options_from_file),options);


if ~options.use_data, data_file = []; end


% ----------------------------------------------------------------------
% Banner

display(sprintf('------------------------------------------------'));
display(sprintf('Parameter balancing for kinetic metabolic models'));
display(sprintf('------------------------------------------------\n'));
display(sprintf('Running parameter balancing\n'))


% ----------------------------------------------------------------------
% Parameter balancing; builds model struct 'network' with balanced kinetic parameters (in field 'kinetics')

pb_options = join_struct(options,struct('parameter_prior_file', prior_file,'parametrisation','all'));

[network, r_mode, r_orig, ~, r_samples, parameter_prior, r_mean, r_std, r_geom_mean, r_geom_std] = parameter_balancing_sbtab(model_file, data_file, pb_options);

% ----------------------------------------------------------------------
% Checks and diagnostic graphics

if options.flag_check,
  parameter_balancing_check(network.kinetics, r_orig, network, parameter_prior,1,1)
end


% ----------------------------------------------------------------------
% convert results (kinetics data structure) to SBtab table struct and save to file

opt                       = struct('write_all_quantities','many','use_sbml_ids',0,'document_name',model_name,'kinetics_mode',r_mode);
opt.more_column_names     = {'UnconstrainedGeometricMean', 'UnconstrainedGeometricStd', 'UnconstrainedMean', 'UnconstrainedStd'};
opt.more_column_data      = {r_geom_mean,r_geom_std,r_mean,r_std};
balanced_parameters_SBtab = modular_rate_law_to_sbtab(network,[],opt);

if length(output_file),
  sbtab_table_save(balanced_parameters_SBtab,struct('filename',output_file));
  display(sprintf('\nWriting output file %s',output_file))
end

% ----------------------------------------------------------------------
% convert posterior samples (kinetics data structure) to SBtab table struct and save to file

if pb_options.n_samples > 0,
  opt                       = struct('write_all_quantities','many','use_sbml_ids',0,'document_name',model_name,'kinetics_mode',r_mode);
  for itt = 1:pb_options.n_samples, opt.more_column_names{itt} = ['Sample', num2str(itt)]; end
  opt.more_column_data      = r_samples;
  sampled_parameter_sets_SBtab = modular_rate_law_to_sbtab(network,[],opt);
  
  if length(output_file),
    output_file_samples = [output_file(1:end-4) '_samples' '.tsv'];
    sbtab_table_save(sampled_parameter_sets_SBtab,struct('filename', output_file_samples));
    display(sprintf('\nWriting sampled parameter sets to output file %s', output_file_samples))
  end
end
