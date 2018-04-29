% -------------------------------------------------------------------------
% Demo script: usage examples for matlab function parameter_balancing_sbtab
% -------------------------------------------------------------------------


% ----------------------------------------------------------------------
% Options

if ~exist('test_example','var'),
  test_example      = 'pfk'; % other choices: 'teusink'; 'hynne'; 'jiang';
  use_pseudo_values = 1; 
  use_data          = 1; 
  flag_check        = 0; 
end

% ----------------------------------------------------------------------
% Load files

model_name = test_example;

[model_file, data_file, prior_file, config_file] = pb_example_files(test_example);

if ~use_data, data_file = []; end

% read options from config file

if 0,
  config     = sbtab_table_load(config_file);
  model_name = sbtab_table_get_element(config,'Value','Option','Model_name');
  use_pseudo_values = strcmp(sbtab_table_get_element(config,'Value','Option','use_pseudos'),'True');
end


% ----------------------------------------------------------------------
% Parameter balancing


display(sprintf('------------------------------------------------'));
display(sprintf('Parameter balancing for metabolic network models'));
display(sprintf('------------------------------------------------\n'));

display('Running parameter balancing')

display(sprintf('  Using model file %s', model_file))

% run parameter balancing; build model struct 'network' with balanced kinetic parameters (in field 'kinetics')

pb_options = struct('parameter_prior_file', prior_file,'use_pseudo_values',use_pseudo_values,'parametrisation','all');

[network, r, r_orig, ~, ~, parameter_prior, r_mean, r_std, r_geom_mean, r_geom_std] = parameter_balancing_sbtab(model_file, data_file, pb_options);

display(sprintf('Parameter balancing finished'));

if flag_check,
  parameter_balancing_check(network.kinetics, r_orig, network, parameter_prior,1,0)
end

% convert kinetics data structure into SBtab table struct

balanced_parameters_SBtab = modular_rate_law_to_sbtab(network,[],struct('write_all_quantities','many','use_sbml_ids',0,'document_name',model_name,'kinetics_mode',r_mean,'kinetics_mean',r_mean,'kinetics_std',r_std,'kinetics_geom_mean',r_geom_mean,'kinetics_geom_std',r_geom_std));

% show table contents (to save the table to disc, insert filename instead of '[]')

balanced_parameters = mytable(sbtab_table_save(balanced_parameters_SBtab),0,[])

%sbtab_table_save(balanced_parameters_SBtab)
