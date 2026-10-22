function  [E_kX_sc,kX_names,kX] = compute_modular_kX_elasticities(kinetic_law,network, c);

% [E_kX_sc,kX_names,kX] = compute_modular_kX_elasticities(kinetic_law,network, c);

eval(default('kinetic_law','cs'));

switch kinetic_law, 
  case 'cs', [E_kX_sc,kX_names,kX] = compute_cs_kX_elasticities(network, c);
  case 'ds', error('compute_modular_kX_elasticities is not implemented for kinetics of type ''ds''');
  case 'ms', [E_kX_sc,kX_names,kX] = compute_ms_kX_elasticities(network, c);
  case 'rp', error('compute_modular_kX_elasticities is not implemented for kinetics of type ''rp''');
  case 'fd', error('compute_modular_kX_elasticities is not implemented for kinetics of type ''fd''');
end
