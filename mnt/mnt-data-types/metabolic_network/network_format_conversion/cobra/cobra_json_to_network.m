function network = cobra_to_network(json_file,options)

% network = cobra_to_network(json_file,options)

options_default = struct('flag_set_external',1);
eval(default('options','struct'));
options = join_struct(options_default,options);

% ----------------------------------

fid = fopen(json_file); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
cobra_model = jsondecode(str);

% ----------------------------------

clear network

for it = 1:length(cobra_model.metabolites),
  network.metabolites{it,1}      = cobra_model.metabolites(it).id;
  network.metabolite_names{it,1} = cobra_model.metabolites(it).name;
  switch cobra_model.metabolites(it).compartment,
    case 'e',
      network.external(it) = 1;
    otherwise
      network.external(it) = 0;
  end
  if isfield(cobra_model.metabolites(it).annotation,'kegg_compound'),
    network.metabolite_KEGGID{it,1} = cobra_model.metabolites(it).annotation.kegg_compound;
  end
  network.metabolite_BIGG{it,1} = cobra_model.metabolites(it).annotation.bigg_metabolite;
end

for it = 1:length(cobra_model.reactions),
  network.actions{it,1}           = cobra_model.reactions{it}.id;
  network.reaction_names{it,1}    = cobra_model.reactions{it}.name;
  if isfield(cobra_model.reactions{it},'annotation'),
  if isfield(cobra_model.reactions{it}.annotation,'kegg_reaction'),
    network.reaction_KEGGID{it,1}    = cobra_model.reactions{it}.annotation.kegg_reaction;
  end
  network.reaction_BIGGID{it,1}    = cobra_model.reactions{it}.annotation.bigg_reaction;
  end
end

network.reversible = ones(size(network.actions));

network.N = zeros(length(network.metabolites), length(network.actions));

for it = 1:length(cobra_model.reactions),
  metlist = fieldnames(cobra_model.reactions{it}.metabolites);
  it
  for it2 = 1:length(metlist),
    my_met = metlist{it2};
    my_stoich = cobra_model.reactions{it}.metabolites.(my_met);
    my_met_ind = label_names(my_met,network.metabolites);
    if my_met_ind==0,
      my_met = my_met(2:end);
      my_met_ind = label_names(my_met,network.metabolites);
    end
    network.N(my_met_ind,it) = my_stoich;
  end
end

network.regulation_matrix = sparse(zeros(size(network.N')));

if options.flag_set_external,
  network = network_set_external(network,1);
end

