function [data, info] = mnt_map_state_data(network, data_files, options)

% [data, info] = mnt_map_state_data(network, data_files, options)
%
% Read metabolic state data (fluxes, metabolite concentrations, enzyme concentrations) from SBtab
% data files and map them to the reactions and metabolites of a model
%
% Inputs
%   network     network structure
%   data_files  SBtab file name, or cell array of file names (the tables can be in the same file or
%               in different files)
%   options     struct with (optional) fields
%     .match_data_by        'ModelElementId' (default), 'KeggId' or 'BiggId': how data rows are matched
%                           to model elements (columns !Reaction / !Compound, or KEGG or BiGG identifiers)
%     .flux_table_id        (default 'FluxData'); [] = no flux data
%     .metabolite_table_id  (default 'MetaboliteConcentrationData'); [] = no concentration data
%     .enzyme_table_id      (default []); e.g. 'EnzymeConcentrationData'
%     .columns_mean         value columns (same names in all tables), e.g. {'4_glucose_Mean'} (default {'Mean'})
%     .columns_std          standard deviation columns (default {}: none)
%     .orient_data          (default 1) convert the data to the reaction orientation of the model first
%                           (sbtab_orient_data_files; the original files are not changed). With 0, the
%                           values are used as they are in the file
%     .model_orientation    (default []) name of the reaction orientation of the model, compared with the
%                           table attribute ReactionOrientation (e.g. 'Model:escherichia-coli-core')
%     .verbose              (default 0)
%
% Outputs
%   data   struct with fields (one row per model reaction or metabolite, one column per value column;
%          NaN for elements without data; empty if the table is not used)
%            .v_data, .v_std   fluxes
%            .c_data, .c_std   metabolite concentrations
%            .e_data, .e_std   enzyme concentrations
%            .v_unit, .c_unit, .e_unit  units (from the tables)
%   info   struct with field .orientation_report (see sbtab_orient_data_files)
%
% Reaction orientation: the data tables may state their reaction orientation (columns
% !ReactionFormula, !ReactionFormula:Identifiers:kegg.compound, table attribute ReactionOrientation);
% see sbtab_orient_data_files. (This handling of reaction orientations will be revised.)
%
% Kinetic and thermodynamic data (e.g. equilibrium constants) are handled by mnt_keq_from_kinetic_data.
% To construct a reference state from these data, see stm_construct_reference_state (STM toolbox).

eval(default('options','struct'));

options_default = struct('match_data_by', 'ModelElementId', 'flux_table_id', 'FluxData', ...
                         'metabolite_table_id', 'MetaboliteConcentrationData', 'enzyme_table_id', [], ...
                         'columns_mean', {{'Mean'}}, 'columns_std', {{}}, ...
                         'orient_data', 1, 'model_orientation', [], 'verbose', 0);
options = join_struct(options_default, options);
if ischar(options.columns_mean), options.columns_mean = {options.columns_mean}; end
if ischar(options.columns_std),  options.columns_std  = {options.columns_std};  end

if ischar(data_files), data_files = {data_files}; end
data_files = unique(data_files(~cellfun(@isempty, data_files)), 'stable');

info.orientation_report = [];
if options.orient_data
  [data_files, info.orientation_report] = sbtab_orient_data_files(network, data_files, ...
       struct('match_data_by', options.match_data_by, 'model_orientation', options.model_orientation, 'verbose', options.verbose));
end

if length(data_files) == 1
  document = sbtab_document_load_from_one(data_files{1});
else
  document = sbtab_document_load(data_files);
end
table_names = sbtab_document_get_table_names(document);

file_options = struct('match_data_by', options.match_data_by, 'columns_mean', {options.columns_mean}, 'columns_std', {options.columns_std});

[data.v_data, data.v_std, data.v_unit] = load_table(network, document, table_names, options.flux_table_id,       'reaction_rate',            file_options);
[data.c_data, data.c_std, data.c_unit] = load_table(network, document, table_names, options.metabolite_table_id, 'metabolite_concentration', file_options);
[data.e_data, data.e_std, data.e_unit] = load_table(network, document, table_names, options.enzyme_table_id,     'enzyme_concentration',     file_options);


% ------------------------------------------------------------

function [x_mean, x_std, unit] = load_table(network, document, table_names, table_id, quantity_type, file_options)

x_mean = []; x_std = []; unit = [];
if isempty(table_id), return; end
if ~any(strcmp(table_id, table_names))
  error('Table %s not found in the data file(s)', table_id);
end
M = load_state_data_for_network(network, sbtab_document_get_table(document, table_id), quantity_type, file_options);
x_mean = M.Mean;
x_std  = M.Std;
unit   = M.Unit;
