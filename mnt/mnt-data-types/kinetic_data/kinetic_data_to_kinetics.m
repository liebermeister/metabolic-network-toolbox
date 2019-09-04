function kinetics = kinetic_data_to_kinetics(kinetics, kinetic_data);

fn = fieldnames(kinetics);
for it = 1:length(fn),
  if isfield(kinetic_data,fn{it}),
    kinetics.(fn{it}) = kinetic_data.(fn{it}).median;
  end
end
