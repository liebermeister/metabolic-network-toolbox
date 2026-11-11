% Conversion of data tables to the reaction orientation of a model
%
% Data tables (SBtab) can state the orientation of their reactions by reaction formulas (columns
% !ReactionFormula and !ReactionFormula:Identifiers:kegg.compound) or by the table attribute
% ReactionOrientation. Before data are mapped to a model, the values are converted to the
% reaction orientation of the model: fluxes and Gibbs free energies change their sign,
% equilibrium constants are inverted, forward and reverse kcat values are exchanged, etc.
% (table quantity_orientation_dependence.tsv in this folder).
%
%   sbtab_orient_data_files                - Convert the data in SBtab files (saved in new files)
%   sbtab_orient_data_table                - Convert one SBtab table
%   sbtab_reaction_orientation             - Orientation of each table row (from formulas or the table attribute)
%   network_reaction_formulas              - Reaction formulas of a network (model compound ids and KEGG ids)
%   sbtab_add_reaction_orientation_columns - Add formulas and the attribute ReactionOrientation to the tables in a data file
%
% Needs the SBtab toolbox with sbtab_table_n_rows and sbtab_table_set_column
