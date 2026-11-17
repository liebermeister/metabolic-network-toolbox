function new_network = network_add_reaction(network,reactants,stoich_coeffs,name,reversible,external)

% NETWORK_ADD_REACTION Add a reaction to a network
%
% new_network = network_add_reaction(network,reactants,stoich_coeffs,name,reversible)
%
% network        network struct
% reactants      metabolite IDs of the reactants; metabolites not yet in the network are added (as external metabolites)
% stoich_coeffs  stoichiometric coefficients
% name           reaction ID (string, or cell array with one string); also used as reaction name
% reversible     (optional, default 1)
%
% The new reaction is appended at the end. All fields of the network are kept (since 2026-11-18; before,
% only the main fields were kept and, for example, metabolite names, compartments and KEGG IDs got lost):
% reaction fields (column vectors of length #reactions, e.g. 'reaction_*', 'enzyme_*', and fields from SBtab
% Reaction tables) are extended by an empty entry ('' or NaN; reaction_names, reaction_NameForPlots, ID: the
% reaction ID); formulae are recomputed.
% Kinetics (network.kinetics) and graphics settings (network.graphics_par) are removed, since they do not
% describe the new reaction.

eval(default('reversible','1'));

if iscell(name), name = name{1}; end

reactants     = column(reactants);
stoich_coeffs = column(stoich_coeffs);

% metabolites that do not yet exist in the network
ll = label_names(reactants,network.metabolites);
additional_metabolites = reactants(find(ll==0));

if length(additional_metabolites)
  additional_external = ones(size(additional_metabolites));
  nn = network_add_metabolites(network,additional_metabolites,additional_external,additional_metabolites);
else
  nn = network;
end

[nr,nm] = network_numbers(nn);

new_network = nn;
ind = label_names(reactants, new_network.metabolites);
new_network.actions                   = [nn.actions; {name}];
new_network.N(:,nr+1)                 = 0;
new_network.N(ind,nr+1)               = stoich_coeffs;
new_network.regulation_matrix(nr+1,:) = 0;
if isfield(nn,'reversible')
  new_network.reversible = [column(nn.reversible); reversible];
else
  new_network.reversible = [ones(nr,1); reversible];
end

% other reaction fields: column vectors of length #reactions (if #metabolites = #reactions: only fields
% 'reaction_*', 'enzyme_*' and fields from SBtab Reaction tables: genes, ID, ReactionFormula, IsReversible)
nm_old = length(network.metabolites);
fn = setdiff(fieldnames(nn), {'metabolites','actions','N','regulation_matrix','reversible','external','formulae'});
for it = 1:length(fn)
  f = fn{it};
  x = nn.(f);
  if size(x,1) ~= nr || size(x,2) ~= 1, continue; end
  is_reaction_field = strncmp(f,'reaction_',9) || strncmp(f,'enzyme_',7) || any(strcmp(f, {'genes','ID','ReactionFormula','IsReversible'}));
  if nm_old == nr && ~is_reaction_field, continue; end
  if iscell(x)
    switch f
      case {'reaction_names','reaction_NameForPlots','ID'}, value = name;
      case 'IsReversible', if reversible, value = 'True'; else, value = 'False'; end
      otherwise, value = '';
    end
    new_network.(f) = [x; {value}];
  else
    new_network.(f) = [x; nan];
  end
end

if isfield(new_network, 'kinetics'),     new_network = rmfield(new_network, 'kinetics');     end
if isfield(new_network, 'graphics_par'), new_network = rmfield(new_network, 'graphics_par'); end

if isfield(nn, 'formulae') || isfield(nn, 'ReactionFormula')
  formulae = network_print_formulae(new_network);
  if isfield(nn, 'formulae'),         new_network.formulae = formulae; end
  if isfield(nn, 'ReactionFormula'),  new_network.ReactionFormula{end} = formulae{end}; end
end
