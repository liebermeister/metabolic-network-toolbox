function [mapping_metabolites, mapping_reactions, covered_metabolites, covered_reactions, shared_metabolites, shared_reactions, network, network_CoHid, network_CoSplit] = embedding_element_mapping(kinetic_models,network,me_options)

% -----------------------------------------------------------------------------------------
% for each kinetic model: 
%  build the mappings mapping_metabolites{it} and mapping_reactions{it} 
%  if elements do not yet exist in the network, create them
%
% Reactions to be mapped must have the same directions!!!
%  (This is checked at the end of this function)

nm_orig = length(network.metabolites);

for it = 1:length(kinetic_models),
  
  %% mapping for metabolites
  mapping_metabolites{it} = label_names(kinetic_models{it}.(me_options.id.metabolites_kinetic_models{it}),network.(me_options.id.metabolites_network));
  ind_nonmapped            = find(mapping_metabolites{it}==0);
  nonmapped_metabolites    = kinetic_models{it}.(me_options.id.metabolites_kinetic_models{it})(ind_nonmapped);
  
  if ind_nonmapped, 
    
    network = network_add_metabolites(network,nonmapped_metabolites,kinetic_models{it}.external(ind_nonmapped));
    display('  Adding metabolites (from kinetic model) to network model:'); 
    fprintf(table(nonmapped_metabolites,1))
    %%display('  Adding them to the network. Metabolites in augmented network:'); 
    %%network.metabolites

    if ~strcmp(me_options.id.metabolites_network,'metabolites'),
      network.(me_options.id.metabolites_network) = [ network.(me_options.id.metabolites_network); kinetic_models{it}.(me_options.id.metabolites_kinetic_models{it})(ind_nonmapped)];
    end
    mapping_metabolites{it} = label_names(kinetic_models{it}.(me_options.id.metabolites_kinetic_models{it}),network.(me_options.id.metabolites_network));  
  end
  
  %% mapping for reactions

  mapping_reactions{it} = label_names(kinetic_models{it}.(me_options.id.reactions_kinetic_models{it}),network.(me_options.id.reactions_network));
  ind_nonmapped          = find(mapping_reactions{it}==0);
  nonmapped_reactions    = kinetic_models{it}.actions(ind_nonmapped);

  if ind_nonmapped,

    nn = network;
    
    for itt = 1:length(ind_nonmapped),
      n = kinetic_models{it}.N(:,ind_nonmapped(itt));
      reactants_kinetic_id = kinetic_models{it}.(me_options.id.metabolites_kinetic_models{it})(find(n));
      ll = label_names(reactants_kinetic_id,nn.(me_options.id.metabolites_network));
      reactants = nn.metabolites(ll);
      stoich_coeffs = n(find(n));
      reversible = kinetic_models{it}.reversible(ind_nonmapped(itt));
      name       = kinetic_models{it}.actions{ind_nonmapped(itt)};
      external   = nn.external(label_names(reactants_kinetic_id,nn.(me_options.id.metabolites_network)));
      nn = network_add_reaction(nn,reactants,stoich_coeffs,name,reversible,external);
      if ~strcmp(me_options.id.reactions_network,'actions'),
        nn.(me_options.id.reactions_network) = [ nn.(me_options.id.reactions_network); kinetic_models{it}.(me_options.id.reactions_kinetic_models{it})(ind_nonmapped(itt))];
      end
    end

    display('  Adding  reactions (from kinetic model) to network model:'); 
    fprintf(table(nonmapped_reactions,1))
    %%display('  Adding them to the network. Reactions in augmented network:'); 
    %%nn.actions

    display('  Resetting field "kinetics" to cs rate law');
    nn.kinetics = set_kinetics(nn,'cs');

    display('  Resetting field "graphics_par"');
    nn = netgraph_make_graph(nn);
    nn = netgraph_read_positions(nn,me_options.position_file,[0,0],1,0);

    mapping_reactions{it} = label_names(kinetic_models{it}.(me_options.id.reactions_kinetic_models{it}),nn.(me_options.id.reactions_network));

    network       = nn;
  end

end

network_CoHid = netgraph_simple_graph(network,me_options.cofactors);
network_CoHid = netgraph_read_positions(network_CoHid,me_options.position_file,[0,0],1,0);
network_CoHid.graphics_par.me_options.omitreactions = me_options.omitreactions;
if isfield(network_CoHid,'metabolite_names'),
  network_CoHid.graphics_par.metnames = escape_underscores(network_CoHid.metabolite_names);
end

network_CoSplit = netgraph_split_metabolites(network,me_options.cofactors);
network_CoSplit = netgraph_read_positions(network_CoSplit,me_options.position_file,[0,0],1,0);
network_CoSplit.graphics_par.me_options.omitreactions = me_options.omitreactions;
if isfield(network_CoSplit,'metabolite_names'),
  network_CoSplit.graphics_par.metnames = escape_underscores(network_CoSplit.metabolite_names);
end


% ---------------------------------------------------------------------
% make vectors showing which network elements are covered by which kinetic models

covered_metabolites = zeros(size(network.metabolites)); 
shared_metabolites  = zeros(size(network.metabolites)); 
covered_reactions   = zeros(size(network.actions)); 
shared_reactions    = zeros(size(network.actions)); 

for  it = length(kinetic_models):-1:1,
  %% bit vectors for newly covered
  d1 = zeros(size(covered_metabolites));
  d2 = zeros(size(covered_reactions));
  d1(mapping_metabolites{it}) = 1;
  d2(mapping_reactions{it})   = 1;
  o1 = d1 .* covered_metabolites;
  o2 = d2 .* covered_reactions;

  covered_metabolites(find(d1)) = it;
  covered_metabolites(find(o1)) = it;
  shared_metabolites(find(o1))  = 1;

  covered_reactions(find(d2))   = it;
  covered_reactions(find(o2))   = it;
  shared_reactions(find(o2))    = 1;

  covered_metabolites(mapping_metabolites{it}) = it;
  covered_reactions(mapping_reactions{it})     = it;
  
end

if find(shared_reactions), 
  display('  Reactions appear in several kinetic models');
  network.actions(find(shared_reactions))
end

if find(shared_metabolites), 
  display('  Metabolites appear in several kinetic models');
  network.metabolites(find(shared_metabolites))
end


% ---------------------------------------------------------------------
% check if reaction stoichiometries match

for it = 1:length(kinetic_models),
  if norm(network.N(mapping_metabolites{it}, mapping_reactions{it}) - kinetic_models{it}.N) > 0,
    error('Conflict between stoichiometries in kinetic models and network');
    else display('  Stoichiometries are consistent')
  end
end


% ---------------------------------------------------------------------
% in principle, metabolites are set internal or external in the combined model
% the same way they were in the network; but all metabolites that are internal 
% in one of the kinetic models are also set internal in the combined model
% and all metabolites that did not exist in the network are set like in the kinetic models

display('  Setting internal metabolites (from kinetic models) internal in network');   
for it = 1:length(kinetic_models),
  network.external(mapping_metabolites{it}(find(kinetic_models{it}.external==0))) = 0;
  ind_new = find(mapping_metabolites{it}>nm_orig);
  network.external(mapping_metabolites{it}(ind_new)) = kinetic_models{it}.external(ind_new);
end
