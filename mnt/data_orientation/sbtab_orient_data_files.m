function [oriented_files, report] = sbtab_orient_data_files(network, data_files, options)

% [oriented_files, report] = sbtab_orient_data_files(network, data_files, options)
%
% Convert the data in SBtab data files to the reaction orientation of a model, and save the
% converted data to new files (in a temporary folder, by default). The original files are not changed.
%
% Each SBtab table with a column !Reaction (or reaction identifier columns) is converted by
% sbtab_orient_data_table (see there). Other tables are copied unchanged.
%
% Inputs
%   network      network structure (see sbtab_reaction_orientation)
%   data_files   file name or cell array of file names (empty entries are kept empty)
%   options      struct with (optional) fields
%     .match_data_by      'ModelElementId' (default), 'KeggId', or 'BiggId'
%     .model_orientation  name of the reaction orientation of the model
%                         (e.g. 'Model:escherichia-coli-core')
%     .output_dir         folder for the converted files (default: [tempdir 'imf-oriented-data'])
%     .verbose            (default 1)
%
% Outputs
%   oriented_files  file names of the converted files (same shape as data_files)
%   report          struct array (one entry per file) with fields .file, .tables (table IDs),
%                   .table_reports (see sbtab_orient_data_table)

eval(default('options','struct'));
options_default = struct('match_data_by', 'ModelElementId', 'model_orientation', [], ...
                         'output_dir', [], 'verbose', 1, 'dependence_file', []);
options = join_struct(options_default, options);

if isempty(options.output_dir)
  options.output_dir = fullfile(tempdir, 'imf-oriented-data');
end
if ~exist(options.output_dir, 'dir'), mkdir(options.output_dir); end

is_char = ischar(data_files);
if is_char, data_files = {data_files}; end

oriented_files = data_files;
report = struct('file', {}, 'tables', {}, 'table_reports', {});

for it_file = 1:length(data_files)

  data_file = data_files{it_file};
  if isempty(data_file), continue; end

  if options.verbose
    display(sprintf('Reaction orientation of data in file %s', data_file));
  end

  my_document = sbtab_document_load_from_one(data_file);
  table_names = fieldnames(my_document.tables);

  report(it_file).file          = data_file;
  report(it_file).tables        = {};
  report(it_file).table_reports = {};

  for it_table = 1:length(table_names)
    my_table = my_document.tables.(table_names{it_table});
    has_reactions = sbtab_table_has_column(my_table, 'Reaction') || ...
        sbtab_table_has_column(my_table, 'Reaction:Identifiers:kegg.reaction') || ...
        sbtab_table_has_column(my_table, 'Reaction:Identifiers:bigg.reaction');
    if ~has_reactions, continue; end
    [my_table, table_report] = sbtab_orient_data_table(network, my_table, options);
    my_document.tables.(table_names{it_table}) = my_table;
    report(it_file).tables{end+1}        = table_names{it_table};
    report(it_file).table_reports{end+1} = table_report;
  end

  % file name of the converted file (with the index of the file, to keep the names unique)
  [~, base_name] = fileparts(data_file);
  oriented_files{it_file} = fullfile(options.output_dir, sprintf('%s_oriented_%d.tsv', base_name, it_file));
  sbtab_document_save_to_one(my_document, oriented_files{it_file}, 0);

end

if is_char, oriented_files = oriented_files{1}; end
