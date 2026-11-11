function [metabolite_bigg_ids, reaction_bigg_ids] = convert_kegg_to_bigg(metabolite_kegg_ids, reaction_kegg_ids, options)

% CONVERT_KEGG_TO_BIGG Convert KEGG ids of metabolites and reactions to BiGG ids
%
% [metabolite_bigg_ids, reaction_bigg_ids] = convert_kegg_to_bigg(metabolite_kegg_ids, reaction_kegg_ids, options)
%
% metabolite_kegg_ids  KEGG compound ids (e.g. 'C00031')
% reaction_kegg_ids    (optional) KEGG reaction ids (e.g. 'R00200')
% options              struct with (optional) fields
%   .id_mapping_file       mapping table (columns as in convert_bigg_to_kegg; default: bigg-kegg-mapping.tsv in this folder)
%   .extra_mapping_files   cell array of SBtab files with additional mappings, e.g. for a specific model;
%                          columns !Identifiers:kegg.compound and !Identifiers:bigg.metabolite (metabolites)
%                          and/or !Identifiers:kegg.reaction and !Identifiers:bigg.reaction (reactions).
%                          Metabolites: these mappings are used first, the main table only for ids not found there.
%                          Reactions: the main table is used first, these mappings only for ids not found there.
%   .metabolites_from_main_table  (default 1) use the main table for metabolites (0: only the extra mappings)
%
% Returned metabolite ids are BiGG ids without compartment (suffixes '_c', '_e', ..., '__L', '__D' removed,
% e.g. 'glc'). Ids without a match are returned as empty strings. The first matching row is used.

if ~exist('reaction_kegg_ids', 'var'), reaction_kegg_ids = {}; end
if ~exist('options', 'var') || ~isstruct(options), options = struct; end
options = join_struct(struct('id_mapping_file', fullfile(fileparts(mfilename('fullpath')), 'bigg-kegg-mapping.tsv'), ...
                             'extra_mapping_files', {{}}, 'metabolites_from_main_table', 1), options);
if ischar(options.extra_mapping_files), options.extra_mapping_files = {options.extra_mapping_files}; end

% extra mappings
extra_kegg_compound = {}; extra_bigg_metabolite = {}; extra_kegg_reaction = {}; extra_bigg_reaction = {};
for it = 1:length(options.extra_mapping_files)
  t = sbtab_table_load(options.extra_mapping_files{it});
  if sbtab_table_has_column(t, 'Identifiers:kegg.compound') && sbtab_table_has_column(t, 'Identifiers:bigg.metabolite')
    extra_kegg_compound   = [extra_kegg_compound;   column(sbtab_table_get_column(t, 'Identifiers:kegg.compound'))];
    extra_bigg_metabolite = [extra_bigg_metabolite; column(sbtab_table_get_column(t, 'Identifiers:bigg.metabolite'))];
  end
  if sbtab_table_has_column(t, 'Identifiers:kegg.reaction') && sbtab_table_has_column(t, 'Identifiers:bigg.reaction')
    extra_kegg_reaction = [extra_kegg_reaction; column(sbtab_table_get_column(t, 'Identifiers:kegg.reaction'))];
    extra_bigg_reaction = [extra_bigg_reaction; column(sbtab_table_get_column(t, 'Identifiers:bigg.reaction'))];
  end
end

% main table
M = bigg_kegg_mapping_load(options.id_mapping_file);

if options.metabolites_from_main_table
  kegg_compound   = [extra_kegg_compound;   M.kegg_compound];
  bigg_metabolite = [extra_bigg_metabolite; M.bigg_metabolite];
else
  kegg_compound   = extra_kegg_compound;
  bigg_metabolite = extra_bigg_metabolite;
end
metabolite_bigg_ids = map_ids(metabolite_kegg_ids, kegg_compound, bigg_kegg_strip_metabolite_ids(bigg_metabolite));
reaction_bigg_ids   = map_ids(reaction_kegg_ids, [M.kegg_reaction; extra_kegg_reaction], [M.bigg_reaction; extra_bigg_reaction]);


function out = map_ids(ids, from, to)

ind  = label_names(ids, from);
out  = repmat({''}, length(ids), 1);
out(ind ~= 0) = to(ind(ind ~= 0));
