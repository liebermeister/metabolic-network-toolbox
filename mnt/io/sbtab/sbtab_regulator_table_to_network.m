function network = sbtab_regulator_table_to_network(network, regulator_table)

% SBTAB_REGULATOR_TABLE_TO_NETWORK Add the regulators of an SBtab Regulator table to a network
%
% network = sbtab_regulator_table_to_network(network, regulator_table)
%
% regulator_table: SBtab table (or filename) of table type Regulator (regulators of
% gene expression, e.g. transcription factors), with the columns
%   !ID                 regulator ID (older files: column !Regulator)
%   !ActivityCondition  (optional) compounds that activate or inhibit the regulator,
%                       as a regulation formula like in !KineticRegulation, e.g. '+ FBP'
%
% New fields of 'network':
%   network.regulators                    regulator IDs (column cell array)
%   network.regulator_activity_condition  entries of column !ActivityCondition (if given)
%   network.regulator_activity_matrix     sparse matrix (regulators x metabolites):
%                                         +1 activation, -1 inhibition (if !ActivityCondition is given)
%
% Compounds that are not in network.metabolites, and compounds with dual or
% unspecified effect, are ignored with a warning.
%
% See also: network_parse_regulation_formula, sbtab_to_network, network_to_sbtab

if ischar(regulator_table)
  regulator_table = sbtab_table_load(regulator_table);
end
regulator_table = sbtab_table_remove_comment_lines(regulator_table);

if sbtab_table_has_column(regulator_table,'ID')
  regulators = sbtab_table_get_column(regulator_table,'ID');
elseif sbtab_table_has_column(regulator_table,'Regulator')
  regulators = sbtab_table_get_column(regulator_table,'Regulator');
else
  warning('Regulator table without column ID; ignoring it');
  return
end
regulators = column(regulators);
network.regulators = regulators;

if ~sbtab_table_has_column(regulator_table,'ActivityCondition')
  return
end

conditions = column(sbtab_table_get_column(regulator_table,'ActivityCondition'));
network.regulator_activity_condition = conditions;

nm = length(network.metabolites);
W  = sparse(length(regulators), nm);
for it = 1:length(regulators)
  [my_signs, my_compounds] = network_parse_regulation_formula(conditions{it});
  if isempty(my_compounds), continue; end
  ind = label_names(my_compounds, network.metabolites);
  if any(ind == 0)
    warning('Regulator %s: unknown compound(s) %s in column ActivityCondition are ignored', regulators{it}, strjoin(my_compounds(ind==0)',', '));
  end
  if any(my_signs(ind ~= 0) == 0)
    warning('Regulator %s ("%s"): compounds with dual or unspecified effect are ignored', regulators{it}, conditions{it});
  end
  ok = (ind ~= 0) & (my_signs(:) ~= 0);
  W(it, ind(ok)) = my_signs(ok);
end
network.regulator_activity_matrix = W;
