% -----------------------------------------
% parameter_balancing_sbtab : usage example
% -----------------------------------------

% set input file names

test_example = 'teusink';

switch test_example,
  case 'pfk',
    model_file  = [ get_pb_directory '/models/pfk.xml'];
    data_file   = [ get_pb_directory '/models/pfk_data.tsv']; 
  case 'teusink',
    model_file  = [ get_pb_directory '/models/teusink.xml'];
    data_file   = [ get_pb_directory '/models/teusink_data.tsv'];    
  case 'hynne',
    model_file  = [ get_pb_directory '/models/hynne.xml'];
    data_file   = [ get_pb_directory '/models/hynne_data.tsv'];    
  case 'jiang',
    model_file  = [ get_pb_directory '/models/jiang.xml'];
    data_file   = [ get_pb_directory '/models/jiang_data.tsv'];
end

prior_file  = [ get_pb_directory '/config/pb_prior_10-2017.tsv'];
config_file = [ get_pb_directory '/config/pb_settings.tsv'];


% set default model name from config file

config     = sbtab_table_load(config_file);

model_name = sbtab_table_get_element(config,'Value','Option','Model_name');

display(sprintf('Balancing the parameters for model "%s"', model_name))

% build model struct 'network' with balanced kinetic parameters (in field 'kinetics')

options = struct('parameter_prior_file', prior_file);

[network, ~, r_orig, ~, ~, parameter_prior, r_mean, r_std] = parameter_balancing_sbtab(model_file, data_file, options);

parameter_balancing_check(network.kinetics, r_orig, network, parameter_prior,1,0)

% convert kinetics data structure into SBtab table struct

balanced_parameters_SBtab = modular_rate_law_to_sbtab(network,[],struct('use_sbml_ids',0,'document_name',model_name,'kinetics_mean',r_mean,'kinetics_std',r_std));

% show table contents (to save the table to disc, insert filename instead of '[]')

mytable(sbtab_table_save(balanced_parameters_SBtab),0,[])

