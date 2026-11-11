function sbtab_add_reaction_orientation_columns(filename, network, model_orientation, options)

% sbtab_add_reaction_orientation_columns(filename, network, model_orientation, options)
%
% Add information about the reaction orientation to the data tables in an SBtab file
% (the file is overwritten). In each table with a column !Reaction (containing the reaction IDs
% of 'network'), this function adds
%   - the columns !ReactionFormula and !ReactionFormula:Identifiers:kegg.compound (formulas
%     of the model reactions, see network_reaction_formulas)
%   - the table attribute ReactionOrientation (value: model_orientation,
%     e.g. 'Model:escherichia-coli-ccm-noor-2016')
%   - optionally, a column !Comment with the text options.comment
% Rows whose reaction is not in the model keep empty formulas.
% With this information, the data can be converted into the reaction orientation of other models
% (see sbtab_orient_data_files).

eval(default('options','struct'));
options = join_struct(struct('comment', []), options);

[formulas, formulas_kegg] = network_reaction_formulas(network);

my_document = sbtab_document_load_from_one(filename);
table_names = fieldnames(my_document.tables);

for it = 1:length(table_names)
  my_table = my_document.tables.(table_names{it});
  if ~sbtab_table_has_column(my_table, 'Reaction'), continue; end
  reaction_ids = column(sbtab_table_get_column(my_table, 'Reaction'));
  [found, ind] = ismember(reaction_ids, network.actions);
  f = repmat({''}, length(reaction_ids), 1); fk = f;
  f(found)  = formulas(ind(found));
  fk(found) = formulas_kegg(ind(found));
  my_table = sbtab_table_set_column(my_table, 'ReactionFormula', f);
  my_table = sbtab_table_set_column(my_table, 'ReactionFormula:Identifiers:kegg.compound', fk);
  if ~isempty(options.comment)
    my_table = sbtab_table_set_column(my_table, 'Comment', repmat({options.comment}, length(reaction_ids), 1));
  end
  my_table.attributes.ReactionOrientation = model_orientation;
  my_document.tables.(table_names{it}) = my_table;
end

sbtab_document_save_to_one(my_document, filename, 0);
