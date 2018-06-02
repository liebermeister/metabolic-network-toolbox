% ---------------------------------------------------------------------------------------------
% Parameter balancing - Example E. coli central carbon metabolism model from Noor et al. (2016)
% Usage example for script 'parameter_balancing_sbtab'
% ---------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------
% Load prepared model and data

model_name = 'E coli central metabolism Noor et al (2016)';

[model_file, data_file] = pb_example_files('ecoli_noor_2016');


% ----------------------------------------------------------------------------
% Set options

pb_options = parameter_balancing_options;

pb_options.v                          = v;
pb_options.enforce_flux_directions    = 0;
pb_options.adjust_to_fluxes           = 0;
pb_options.preferred_data_element_ids = 'sbml';

% ----------------------------------------------------------------------------
% Balance the model parameters

[network, r, r_orig, ~, ~, parameter_prior] = parameter_balancing_sbtab(model_file, data_file, pb_options);

%----------------------------------------------------- 
% An alternative (using an SBtab model file) would be:
% [network, v, c_data, u_data, conc_min, conc_max, positions, warnings] = load_model_and_data_sbtab(sbtab_model_file, data_file, pb_options);
% [r, r_orig, ~, ~, parameter_prior] = parameter_balancing_kinetic(network, data_file, pb_options);
%----------------------------------------------------- 

parameter_balancing_check(r, r_orig, network, parameter_prior,1,1)


% ----------------------------------------------------------------------------
% Convert "kinetics" data structure (used in metabolic network toolbox models)
% into sbtab data structure

network.kinetics = r;

balanced_parameters_sbtab = modular_rate_law_to_sbtab(network,[],struct('use_sbml_ids',0,'document_name',model_name,'write_concentrations',1,'write_enzyme_concentrations',1));


% ----------------------------------------------------------------------------
% Show table contents (to save the table to disc, replace '[]' by output filename)

%mytable(sbtab_table_save(balanced_parameters_sbtab),0,[])

