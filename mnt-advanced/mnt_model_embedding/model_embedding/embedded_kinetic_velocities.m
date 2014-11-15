function v = embedded_kinetic_velocities(c,network_combined,kinetics)

eval(default('kinetics','network_combined.kinetics'));

% evaluate rates for the network model

v = network_velocities(c,network_combined,kinetics.kinetics_network);

% evaluate rates for the kinetic models and insert them
% if there are overlaps between kinetic models, the first one of them wins

for it = length(network_combined.kinetics.kinetic_models):-1:1,
  kk = network_combined.kinetics.kinetic_models{it};
  my_c = c(network_combined.kinetics.mapping_metabolites{it});
  vv = network_velocities(my_c,kk);
  if isfield(kinetics,'enzyme_adjustment'),
    vv =  vv .* kinetics.enzyme_adjustment{it};
  end
  v(network_combined.kinetics.mapping_reactions{it}) = vv;
end

