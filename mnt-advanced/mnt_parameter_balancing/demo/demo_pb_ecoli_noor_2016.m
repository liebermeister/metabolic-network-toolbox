% ---------------------------------------------------------------------------------------------
% Parameter balancing - Example E. coli central carbon metabolism model from Noor et al. (2016)
% ---------------------------------------------------------------------------------------------

model_name  = 'E coli central metabolism Noor et al (2016)';
model_file  = [ get_pb_directory '/models/ecoli_noor_2016.tsv'];
data_file   = [ get_pb_directory '/models/ecoli_noor_2016_data.tsv']; 

[network, v, c_data, u_data, conc_min, conc_max, positions, warnings] = load_model_and_data_sbtab(model_file,data_file);

options = parameter_balancing_default_options;

options.enforce_flux_directions       = 1;
options.postprocessing_enforce_fluxes = 1;

% compute balanced model parameters

[r, r_orig, ~, ~, parameter_prior] = parameter_balancing_kinetic_metabolic(network, data_file, options, v);

parameter_balancing_check(r, r_orig, network, parameter_prior,1,1)

% convert "kinetics" data structure (used in metabolic network toolbox models) into sbtab data structure

network.kinetics = r;

balanced_parameters_sbtab = modular_rate_law_to_sbtab(network,[],struct('use_sbml_ids',0,'document_name',model_name,'write_concentrations',1,'write_enzyme_concentrations',1));

% show table contents 
% (to save the table to disc, replace '[]' by filename of output file)

mytable(sbtab_table_save(balanced_parameters_sbtab),0,[])

