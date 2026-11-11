function network = network_add_kinetic_regulation(network, regulation, options)

% NETWORK_ADD_KINETIC_REGULATION Add kinetic regulation arrows (from an SBtab Reaction table) to a network
%
% network = network_add_kinetic_regulation(network, regulation, options)
%
% The regulation is given as an SBtab Reaction table with the columns !ID (reaction ids of the
% network) and !KineticRegulation (regulation formula, e.g. '- NADPH' or '+ ADP - PEP'; see
% network_parse_regulation_formula). Other columns (e.g. !ReactionFormula) are not needed and are ignored.
% Example file:
%
%   !!!SBtab Document='Regulation'
%   !!SBtab TableID='Reaction' TableType='Reaction'
%   !ID           !KineticRegulation
%   ZWF_R00835    - NADPH
%
% Inputs
%   network     network structure (field regulation_matrix: reactions x metabolites)
%   regulation  SBtab file name, SBtab table object, or SBtab document object (table with TableID
%               'Reaction', or else the first table of TableType 'Reaction')
%   options     struct with (optional) fields
%     .replace  (default 0) 0: the regulation arrows are added to the existing ones;
%               1: the existing regulation of the listed reactions is replaced
%
% Rows with an empty regulation formula are ignored (with replace = 1, they remove the regulation of
% the reaction). Regulators with dual or unknown effect are ignored with a warning. Unknown reaction
% or compound ids give an error.

eval(default('options', 'struct'));
options = join_struct(struct('replace', 0), options);

table = sbtab_table_from_input(regulation, {'Reaction'}, 'Reaction');
if ~sbtab_table_has_column(table, 'ID') || ~sbtab_table_has_column(table, 'KineticRegulation')
  error('The regulation table needs the columns !ID and !KineticRegulation');
end
reaction_ids = sbtab_table_get_column(table, 'ID');
formulas     = sbtab_table_get_column(table, 'KineticRegulation');

[nm, nr] = size(network.N);
if ~isfield(network, 'regulation_matrix') || isempty(network.regulation_matrix)
  network.regulation_matrix = zeros(nr, nm);
end

ind_reactions = label_names(reaction_ids, network.actions);
if any(ind_reactions == 0)
  error('Unknown reaction(s) in the regulation table: %s', strjoin(column(reaction_ids(ind_reactions == 0))', ', '));
end

for it = 1:length(reaction_ids)
  if options.replace
    network.regulation_matrix(ind_reactions(it), :) = 0;
  end
  [signs, compounds] = network_parse_regulation_formula(formulas{it});
  if isempty(compounds), continue; end
  ind_compounds = label_names(compounds, network.metabolites);
  if any(ind_compounds == 0)
    error('Unknown compound(s) %s in the regulation of reaction %s', strjoin(column(compounds(ind_compounds == 0))', ', '), reaction_ids{it});
  end
  if any(signs == 0)
    warning('Reaction %s ("%s"): regulators with dual or unknown effect are ignored', reaction_ids{it}, formulas{it});
  end
  network.regulation_matrix(ind_reactions(it), ind_compounds(signs ~= 0)) = signs(signs ~= 0);
end
