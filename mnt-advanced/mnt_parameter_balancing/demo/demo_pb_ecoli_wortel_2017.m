% -----------------------------------------------------------------------------------------
% Parameter balancing applied to E coli central carbon metabolism model (Wortel et al 2017)
%
% Note that the results differ from the model parameters used in the article
% (here, standard settings are used; in the article, the priors were modified)
% -----------------------------------------------------------------------------------------

model_name  = ' E coli central carbon metabolism Wortel et al (2017)';
model_file  = [ get_pb_directory '/models/ecoli_wortel_2017.xml'];
data_file   = [ get_pb_directory '/models/ecoli_wortel_2017_data.tsv'];
prior_file  = [ get_pb_directory '/config/pb_prior_10-2017.tsv'];


display(sprintf('Balancing the parameters for model "%s"', model_name))

% build model struct 'network' with balanced kinetic parameters (in field 'kinetics')

options = struct('parameter_prior_file', prior_file);

[network, r, r_orig, ~, ~, parameter_prior] = parameter_balancing_sbtab(model_file, data_file, options);

parameter_balancing_check(r, r_orig, network, parameter_prior)

% convert kinetics data structure into SBtab table struct

balanced_parameters_SBtab = modular_rate_law_to_sbtab(network,[],struct('use_sbml_ids',0,'document_name',model_name));

% show table contents (to save the table to disc, insert filename instead of '[]')

mytable(sbtab_table_save(balanced_parameters_SBtab),0,[])

