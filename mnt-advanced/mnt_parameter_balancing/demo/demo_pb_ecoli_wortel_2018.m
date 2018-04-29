% -----------------------------------------------------------------------------------------
% Parameter balancing applied to E coli central carbon metabolism model (Wortel et al 2017)
% -----------------------------------------------------------------------------------------

model_name    = ' E coli central carbon metabolism Wortel et al (2017)';
model_file    = [ get_pb_directory '/models/ecoli_wortel_2017.xml'];
data_file     = [ get_pb_directory '/models/ecoli_wortel_2017_data.tsv'];
pb_prior_file = [ get_pb_directory '/config/pb_prior.tsv'];


display(sprintf('Balancing the parameters for model "%s"', model_name))

% build model struct 'network' with balanced kinetic parameters (in field 'kinetics')

pb_options                      = parameter_balancing_default_options;
pb_options.parameter_prior_file = prior_file;
pb_options.Keq_upper            = 10^10;
pb_options.kcat_prior_median    = 100; % 1/s;  default=10
pb_options.kcat_prior_log10_std = 1.5; % 1/s;  default = 0.2; vorher mal 0.0002
pb_options.parameter_prior_file = pb_prior_file;
pb_options.GFE_fixed            = 0;
pb_options.use_pseudo_values    = 1;

[network, r, r_orig, ~, ~, parameter_prior] = parameter_balancing_sbtab(model_file, data_file, pb_options);

parameter_balancing_check(r, r_orig, network, parameter_prior)

% convert kinetics data structure into SBtab table struct

balanced_parameters_SBtab = modular_rate_law_to_sbtab(network,[],struct('use_sbml_ids',0,'document_name',model_name));

% show table contents (to save the table to disc, insert filename instead of '[]')

mytable(sbtab_table_save(balanced_parameters_SBtab),0,[])

