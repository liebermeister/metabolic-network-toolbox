function pb_options = parameter_balancing_update_options(pb_options)
  
display('Updating the parameters for parameter_balancing');

pb_options_default = parameter_balancing_default_options;
pb_options         = join_struct(pb_options_default,pb_options);
