function [network,c,v,dmu,keq,theta] = sbtab_network_state_load(filename)

% [network,c,v,dmu,keq,theta] = sbtab_network_state_load(filename)
%
% Load 'network' data structure, concentration vector 'c', and flux vector 'v' from SBtab file
% (and, optionally, reaction Gibbs free energies 'dmu' and equilibrium constants 'keq'; empty if not in the file)
%
% filename: SBtab file name or SBtab document object
% Tables: MetaboliteConcentration, Flux, ReactionGibbsFreeEnergy, EquilibriumConstant (TableIDs).
% If the table MetaboliteConcentration or Flux is missing, c or v is a vector of NaN values (with a warning).

sbtab = sbtab_get(filename);

network = sbtab_to_network(sbtab);
table_names = sbtab_document_get_table_names(sbtab);

c = nan(length(network.metabolites),1);
if any(strcmp('MetaboliteConcentration', table_names))
  sbtab_c = sbtab_to_struct(sbtab_document_get_table(sbtab,'MetaboliteConcentration'),'column');
  ll = label_names(network.metabolites,sbtab_c.Metabolite);
  c(ll~=0) = cell_string2num(sbtab_c.Value(ll(ll~=0)));
else
  warning('No metabolite concentrations found in file');
end

v = nan(length(network.actions),1);
if any(strcmp('Flux', table_names))
  sbtab_v = sbtab_to_struct(sbtab_document_get_table(sbtab,'Flux'), 'column');
  ll = label_names(network.actions,sbtab_v.Reaction);
  v(ll~=0) = cell_string2num(sbtab_v.Value(ll(ll~=0)));
else
  warning('No fluxes found in file');
end

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

theta = [];
if nargout > 5 && ~isempty(keq)
  theta = log(keq) - network.N' * log(c);
end
