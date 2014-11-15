function nn = model_embedding_kinetics_combined2flat(network_combined)

kk = network_combined.kinetics;

% kinetics type must be 'embedded_kinetic_models'

if ~strcmp(kk.type, 'embedded_kinetic_models'),
  error('Wrong kinetics type');
end

% all metabolite names must match

for it = 1:length(kk.mapping_metabolites),
   map = kk.mapping_metabolites{it};
   met_kinetic_model = kk.kinetic_models{it}.metabolites;
   met_network_model = network_combined.metabolites(map);
   if ~prod(double(strcmp(met_kinetic_model, met_network_model))),
     error('Metabolite IDs do not match between network model and embedded kinetic model');
   end
end

% convert network kinetics to type 'kinetic_strings'

kk.kinetics_network = kinetics_convert_to_strings(network_combined,kk.kinetics_network);

% convert kinetic submodel kinetics to type 'kinetic_strings'

for it = 1:length(kk.kinetic_models),
   kk.kinetic_models{it}.kinetics = kinetics_convert_to_strings(kk.kinetic_models{it}, kk.kinetic_models{it}.kinetics);
  if length(kk.kinetic_models{it}.kinetics.parameters),
    error('Cannot convert kinetics with global parameters');
  end
end

% make new kinetics field for entire model
% first insert 'kinetic_strings' type kinetics for entire network ..

kkk.type             = 'kinetic_strings';
kkk.reactions        = kk.kinetics_network.reactions;
kkk.parameters       = {};
kkk.parameter_values = [];

% .. then insert rate laws from submodels

for it = 1:length(kk.kinetic_models),
  ind_reactions = kk.mapping_reactions{it};
  my_rate_laws  = kk.kinetic_models{it}.kinetics.reactions;
  my_enzyme_adjustments = kk.enzyme_adjustment{1};
  %% don't forget to adjust enzyme levels!!
  for itt = 1:length( my_rate_laws),
   my_rate_laws{itt}.string = sprintf('EnzymeAdjustment * ( %s )', my_rate_laws{itt}.string);
   my_rate_laws{itt}.parameters = [my_rate_laws{itt}.parameters; 'EnzymeAdjustment'];
   my_rate_laws{itt}.parameter_values = [my_rate_laws{itt}.parameter_values; my_enzyme_adjustments(itt)];
  end
  kkk.reactions(ind_reactions) = my_rate_laws;
end

nn = network_combined;
nn.kinetics = kkk;