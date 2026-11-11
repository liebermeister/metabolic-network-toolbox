function [v,network,ind_keep] = network_remove_zero_fluxes(v,network)

% [v,network,ind_keep] = network_remove_zero_fluxes(v,network)
%
% Remove reactions with zero flux from a flux vector and a network struct (data structure from MNT)
% (metabolites are kept)
% ind_keep: logical vector marking the reactions that were kept
%
% Kinetics (network.kinetics) are reduced to the remaining reactions for the kinetics types
% 'cs', 'ms', 'rp', 'ma', 'fm' (modular rate laws), 'mass-action', and 'kinetic_strings'.
% In other cases (other kinetics types, no field 'type', or field sizes that do not match the
% network), the kinetics are removed, with a warning.

[nm, nr] = size(network.N);

ind_keep = v ~= 0;
indr     = find(ind_keep);
v = v(ind_keep);

network.N = network.N(:,ind_keep);
if isfield(network,'regulation_matrix')
  network.regulation_matrix = network.regulation_matrix(ind_keep,:);
end
network.reversible = network.reversible(ind_keep);
network.actions    = network.actions(ind_keep);

reaction_fields = {'reaction_NameForPlots', 'reaction_names', 'reaction_BIGGID', 'reaction_KEGGID', 'genes', 'formulae', ...
                   'ID', 'ReactionFormula', 'IsReversible', 'enzyme_mass', 'EC'};
for it = 1:length(reaction_fields)
  fn = reaction_fields{it};
  if isfield(network,fn) && numel(network.(fn)) == nr
    network.(fn) = network.(fn)(ind_keep);
  end
end

% kinetics

if isfield(network,'kinetics')
  [kinetics, ok] = reduce_kinetics(network.kinetics, nm, nr, indr);
  if ok
    network.kinetics = kinetics;
  else
    warning('network_remove_zero_fluxes: network.kinetics could not be reduced to the remaining reactions (unsupported kinetics type, no field ''type'', or sizes that do not match the network) and was removed');
    network = rmfield(network,'kinetics');
  end
end


% ------------------------------------------------------------------

function [kinetics, ok] = reduce_kinetics(kinetics, nm, nr, indr)

ok = 0;
if ~isfield(kinetics,'type'), return; end

try
  switch kinetics.type
    case {'cs','ms','rp','ma','fm'}
      for fn = {'u','KV','Keq'}
        if isfield(kinetics,fn{1}) && numel(kinetics.(fn{1})) ~= nr, return; end
      end
      for fn = {'KM','KA','KI'}
        if isfield(kinetics,fn{1}) && ~isequal(size(kinetics.(fn{1})), [nr, nm]), return; end
      end
      kinetics = modular_reduce_to_subnetwork(kinetics, (1:nm)', indr);
    case 'mass-action'
      if numel(kinetics.k_fwd) ~= nr || numel(kinetics.k_bwd) ~= nr, return; end
      kinetics.k_fwd = kinetics.k_fwd(indr);
      kinetics.k_bwd = kinetics.k_bwd(indr);
    case 'kinetic_strings'
      if numel(kinetics.reactions) ~= nr, return; end
      kinetics = submodel_kinetic_strings(kinetics, nm, (1:nm)', indr);
    otherwise
      return
  end
  ok = 1;
catch
  ok = 0;
end
