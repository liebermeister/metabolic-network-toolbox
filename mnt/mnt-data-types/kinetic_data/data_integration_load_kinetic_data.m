function kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, data_file, use_sbml_ids, use_kegg_ids, flag_invent_std, verbose)

% wrapper function, used for backward compatibility
  
options.use_sbml_ids    = use_sbml_ids;
options.use_kegg_ids    = use_kegg_ids;
options.flag_invent_std = flag_invent_std;
options.verbose         = verbose

kinetic_data = kinetic_data_load(data_quantities, parameter_prior, network, data_file, options)