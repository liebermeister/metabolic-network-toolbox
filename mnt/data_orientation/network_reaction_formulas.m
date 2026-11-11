function [formulas, formulas_kegg] = network_reaction_formulas(network)

% [formulas, formulas_kegg] = network_reaction_formulas(network)
%
% Reaction formulas of a network model, written with the compound IDs of the model
% (network.metabolites) and with KEGG compound IDs (network.metabolite_KEGGID; compounds
% without KEGG ID are written with their model ID)
%
% Example: 'ADP + _3_Phospho_D_glyceroyl_phosphate <=> ATP + _3_Phospho_D_glycerate'
%          'C00008 + C00236 <=> C00002 + C00197'

kegg_ids = network.metabolites;
if isfield(network, 'metabolite_KEGGID')
  for it = 1:length(kegg_ids)
    if ~isempty(network.metabolite_KEGGID{it}), kegg_ids{it} = network.metabolite_KEGGID{it}; end
  end
end

n_reactions   = size(network.N, 2);
formulas      = cell(n_reactions, 1);
formulas_kegg = cell(n_reactions, 1);

for j = 1:n_reactions
  formulas{j}      = one_formula(network.N(:,j), network.metabolites);
  formulas_kegg{j} = one_formula(network.N(:,j), kegg_ids);
end


% ------------------------------------------------------------------

function f = one_formula(n, ids)

sides = {};
for s = [-1, 1]
  ind   = find(sign(n) == s);
  terms = {};
  for i = column(ind)'
    c = abs(full(n(i)));
    if c == 1
      terms{end+1} = ids{i};
    else
      terms{end+1} = sprintf('%g %s', c, ids{i});
    end
  end
  sides{end+1} = strjoin(terms, ' + ');
end
f = [sides{1} ' <=> ' sides{2}];
