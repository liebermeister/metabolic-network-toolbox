% -----------------------------------------
% parameter_balancing_sbtab : usage example
% -----------------------------------------

% set input file names

test_example = 'pfk';

switch test_example,
  case 'pfk',
    model_file  = [ get_pb_directory '/demo/models/pfk.xml'];
    data_file   = [ get_pb_directory '/demo/models/pfk_parameters.tsv']; 
  case 'teusink',
    model_file  = [ get_pb_directory '/demo/models/teusink.xml'];
    data_file   = [ get_pb_directory '/demo/models/teusink_parameters.tsv'];    
  case 'hynne',
    model_file  = [ get_pb_directory '/demo/models/hynne.xml'];
    data_file   = [ get_pb_directory '/demo/models/hynne_parameters.tsv'];    
  case 'jiang',
    model_file  = [ get_pb_directory '/demo/models/jiang.xml'];
    data_file   = [ get_pb_directory '/demo/models/jiang_parameters.tsv'];
end


prior_file  = [ get_pb_directory '/demo/config/prior.tsv'];
config_file = [ get_pb_directory '/demo/config/config.tsv'];


% set default model name from config file

config     = sbtab_table_load(config_file);

model_name = sbtab_table_get_element(config,'Value','Option','Model_name');

display(sprintf('Balancing the parameters for model "%s"', model_name))

% build model struct 'network' with balanced kinetic parameters (in field 'kinetics')

options = struct('parameter_prior_file', prior_file);

[network, ~, ~, ~, ~, ~, r_mean, r_std] = parameter_balancing_sbtab(model_file, data_file, options);


% convert kinetics data structure into SBtab table struct

balanced_parameters_SBtab = modular_rate_law_to_sbtab(network,[],struct('use_sbml_ids',0,'document_name',model_name,'kinetics_mean',r_mean,'kinetics_std',r_std));


% show table contents (to save the table to disc, insert filename instead of '[]')

mytable(sbtab_table_save(balanced_parameters_SBtab),0,[])

