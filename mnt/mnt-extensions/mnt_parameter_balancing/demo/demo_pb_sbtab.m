% --------------------------------------------------
% Parameter balancing demo script
% Usage example for matlab function 'parameter_balancing'
%
% Example models: example models the models/ subdirectory
%  'pfk'                - example from Parameter balancing website
%  'teusink'            - example from Parameter balancing website
%  'hynne'              - example from Parameter balancing website
%  'hynne_inh'          - example from Parameter balancing website
%  'jiang'              - example from Parameter balancing website
%  'ecoli_noor_2016'    - E. coli model, Noor et al 2016
%  'ecoli_wortel_2018'  - E. coli model, Wortel et al 2018
%
% (Filenames are provided by the function pb_example_files)
% --------------------------------------------------

% ----------------------------------------------------------------------
% Choose model from "models" directory

test_example = 'ecoli_noor_2016';
%test_example = 'ecoli_wortel_2018';
%test_example = 'pfk';
%test_example = 'teusink';
%test_example = 'hynne';
%test_example = 'hynne_inh';
%test_example = 'jiang';

% ----------------------------------------------------------------------
% Options

eval(default('options','struct'));

options.preferred_data_element_ids = 'sbml';
options.flag_check                 = 1; 

% Further options, overriding options in options file (optional)
% options.use_pseudo_values = 1; 
% options.use_data          = 1; 
% options.n_samples         = 0; 

% ----------------------------------------------------------------------
% Load model and data files

[model_name, model_file, data_file, prior_file, options_file] = pb_example_files(test_example);

% To save result files, replace '[]' by [PATHNAME]/[FILENAME]
output_file = [];


% ----------------------------------------------------------------------
% Run parameter balancing

close all

balanced_parameters_SBtab = parameter_balancing(model_file, output_file, data_file, prior_file, options_file, model_name, options);

%Show results (balanced parameters) as a table:
% mytable(sbtab_table_save(balanced_parameters_SBtab),0,[])

