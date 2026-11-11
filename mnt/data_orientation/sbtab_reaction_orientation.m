function [orientation, ind_reaction, method] = sbtab_reaction_orientation(network, data_table, options)

% [orientation, ind_reaction, method] = sbtab_reaction_orientation(network, data_table, options)
%
% Determine, for each row of an SBtab data table, whether the reaction is oriented in the same
% way as in the model, or in the opposite way.
%
% Inputs
%   network      network structure (fields actions, metabolites, N; optional: reaction_KEGGID,
%                reaction_BIGGID, metabolite_KEGGID)
%   data_table   SBtab table object with a column !Reaction (or reaction identifier columns)
%   options      struct with (optional) fields
%     .match_data_by      'ModelElementId' (default), 'KeggId', or 'BiggId': how table rows are
%                         matched to model reactions (as in sbtab_load_network_model_and_data)
%     .model_orientation  name of the reaction orientation of the model (e.g.
%                         'Model:escherichia-coli-ccm-noor-2016'); compared with the table
%                         attribute ReactionOrientation
%
% Outputs (one entry per table row)
%   orientation   +1 (same orientation as in the model), -1 (opposite orientation), or NaN (unknown,
%                 or the row does not refer to a model reaction)
%   ind_reaction  index of the model reaction (NaN if the row does not match a model reaction)
%   method        how the orientation was determined: 'formula-kegg', 'formula', 'table-attribute',
%                 or '' (unknown)
%
% The orientation of a row is determined (in this order of priority) from
%   1. the column !ReactionFormula:Identifiers:kegg.compound (reaction formula with KEGG compound
%      IDs), compared with the model reaction via network.metabolite_KEGGID
%   2. the column !ReactionFormula (reaction formula with the compound IDs of the model)
%   3. the table attribute ReactionOrientation: if it equals options.model_orientation, all rows
%      have the orientation of the model
% Formulas are compared by the compounds that occur in both formulas (H+ and H2O are ignored):
% compounds on the same side vote for the same orientation, compounds on opposite sides vote
% for the opposite orientation.

eval(default('options','struct'));
options_default = struct('match_data_by', 'ModelElementId', 'model_orientation', []);
options = join_struct(options_default, options);

n_rows = sbtab_table_n_rows(data_table);
orientation  = nan(n_rows,1);
ind_reaction = nan(n_rows,1);
method       = repmat({''}, n_rows, 1);

% -------------------------------------
% match rows to model reactions

switch options.match_data_by
  case 'ModelElementId',
    data_ids  = get_column(data_table, 'Reaction');
    model_ids = network.actions;
  case 'KeggId',
    data_ids  = get_column(data_table, 'Reaction:Identifiers:kegg.reaction');
    model_ids = network.reaction_KEGGID;
  case 'BiggId',
    data_ids  = get_column(data_table, 'Reaction:Identifiers:bigg.reaction');
    model_ids = network.reaction_BIGGID;
  otherwise,
    error('Unknown option match_data_by: %s', options.match_data_by);
end

if isempty(data_ids), return; end

for it = 1:n_rows
  if ~isempty(data_ids{it})
    ll = find(strcmp(data_ids{it}, model_ids));
    if length(ll) == 1, ind_reaction(it) = ll; end
  end
end

% -------------------------------------
% orientation from formulas or from the table attribute

formula_kegg = get_column(data_table, 'ReactionFormula:Identifiers:kegg.compound');
formula      = get_column(data_table, 'ReactionFormula');

table_orientation = '';
if isfield(data_table.attributes, 'ReactionOrientation')
  table_orientation = data_table.attributes.ReactionOrientation;
end

has_kegg = isfield(network, 'metabolite_KEGGID');
ignore_kegg = {'C00080', 'C00001', 'C18714'}; % H+, H2O
ignore_ids  = {'H2O', 'h2o', 'H', 'h', 'Hplus', 'h2o[c]', 'h2o[e]', 'h[c]', 'h[e]'};

for it = find(isfinite(ind_reaction))'
  j = ind_reaction(it);
  o = nan;
  if has_kegg && ~isempty(formula_kegg) && ~isempty(formula_kegg{it})
    o = compare_formula(formula_kegg{it}, network.N(:,j), network.metabolite_KEGGID, ignore_kegg);
    if isfinite(o), method{it} = 'formula-kegg'; end
  end
  if isnan(o) && ~isempty(formula) && ~isempty(formula{it})
    o = compare_formula(formula{it}, network.N(:,j), network.metabolites, ignore_ids);
    if isfinite(o), method{it} = 'formula'; end
  end
  if isnan(o) && ~isempty(table_orientation) && ~isempty(options.model_orientation) ...
        && strcmp(table_orientation, options.model_orientation)
    o = 1; method{it} = 'table-attribute';
  end
  orientation(it) = o;
end


% ------------------------------------------------------------------

function o = compare_formula(formula_string, n_column, model_compound_ids, ignore_ids)

% compare a reaction formula (string) with a column of the stoichiometric matrix
% o = 1 (same orientation), -1 (opposite orientation), nan (no decision)

o = nan;
sides = regexp(formula_string, '<=>|<->|->|=', 'split');
if length(sides) ~= 2, return; end

score = 0;
side_signs = [-1, 1];
for is = 1:2
  terms = strtrim(strsplit(sides{is}, ' + '));
  for t = 1:length(terms)
    term = strtrim(terms{t});
    if isempty(term), continue; end
    parts = strsplit(term, ' ');
    compound = parts{end};
    if sum(strcmp(compound, ignore_ids)), continue; end
    ll = find(strcmp(compound, model_compound_ids));
    if isempty(ll), continue; end
    model_stoich = sum(n_column(ll));
    if model_stoich ~= 0
      score = score + side_signs(is) * sign(model_stoich);
    end
  end
end

if score > 0, o = 1; elseif score < 0, o = -1; end


% ------------------------------------------------------------------

function c = get_column(data_table, column_name)

c = sbtab_table_get_column(data_table, column_name, 0, 0);
if ~isempty(c), c = column(c); end
