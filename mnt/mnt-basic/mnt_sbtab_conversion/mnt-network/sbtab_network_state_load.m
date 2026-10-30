function [network,c,v,dmu,keq] = sbtab_network_state_load(filename)

% [network,c,v,dmu,keq] = sbtab_network_state_load(filename)
%
% Load 'network' data structure, concentration vector 'c', and flux vector 'v' from SBtab file
% (and, optionally, reaction Gibbs free energies 'dmu' and equilibrium constants 'keq'; empty if not in the file)

sbtab = sbtab_document_load(filename);

network = sbtab_to_network(sbtab);
sbtab_c = sbtab_to_struct(sbtab_document_get_table(sbtab,'MetaboliteConcentration'),'column');
sbtab_v = sbtab_to_struct(sbtab_document_get_table(sbtab,'Flux'), 'column');

ll = label_names(network.metabolites,sbtab_c.Metabolite);
c  = cell_string2num(sbtab_c.Value(ll));

ll = label_names(network.actions,sbtab_v.Reaction);
v  = cell_string2num(sbtab_v.Value(ll));

dmu = [];
keq = [];

if nargout > 3
  if find(strcmp('ReactionGibbsFreeEnergy',sbtab_document_get_table_names(sbtab)))
    sbtab_dmu = sbtab_to_struct(sbtab_document_get_table(sbtab,'ReactionGibbsFreeEnergy'),'column');
    ll = label_names(network.actions,sbtab_dmu.Reaction);
    dmu  = cell_string2num(sbtab_dmu.Value(ll));
  else
    warning('No Reaction Gibbs free energy values found in file');
  end
end

if nargout > 4
  if find(strcmp('EquilibriumConstant',sbtab_document_get_table_names(sbtab)))
    sbtab_dmu = sbtab_to_struct(sbtab_document_get_table(sbtab,'EquilibriumConstant'),'column');
    ll = label_names(network.actions,sbtab_dmu.Reaction);
    keq  = cell_string2num(sbtab_dmu.Value(ll));
  else
    warning('No equilibrium constant values found in file');
  end
end
