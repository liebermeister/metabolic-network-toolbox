function [metabolite_kegg_ids, reaction_kegg_ids] = convert_bigg_to_kegg(metabolite_bigg_ids, reaction_bigg_ids, id_mapping_file)

% CONVERT_BIGG_TO_KEGG Convert BiGG ids of metabolites and reactions to KEGG ids
%
% [metabolite_kegg_ids, reaction_kegg_ids] = convert_bigg_to_kegg(metabolite_bigg_ids, reaction_bigg_ids, id_mapping_file)
%
% metabolite_bigg_ids  BiGG metabolite ids without compartment suffix, e.g. 'glc' or 'glc[e]' ('[c]' and '[e]'
%                      are ignored); they are compared with the ids of the table without the suffixes
%                      '_c', '_e', ..., '__L', '__D' (see bigg_kegg_strip_metabolite_ids)
% reaction_bigg_ids    (optional) BiGG reaction ids
% id_mapping_file      (optional) mapping table with the columns !Identifiers:bigg.reaction,
%                      !Identifiers:bigg.metabolite, !Identifiers:kegg.reaction, !Identifiers:kegg.compound
%                      (default: bigg-kegg-mapping.tsv in this folder)
%
% Ids without a match are returned as empty strings. The first matching row of the table is used (also if
% its KEGG id is empty).

if ~exist('reaction_bigg_ids', 'var'), reaction_bigg_ids = {}; end
if ~exist('id_mapping_file', 'var') || isempty(id_mapping_file)
  id_mapping_file = fullfile(fileparts(mfilename('fullpath')), 'bigg-kegg-mapping.tsv');
end

M = bigg_kegg_mapping_load(id_mapping_file);

metabolite_bigg_ids = strrep(strrep(metabolite_bigg_ids, '[c]', ''), '[e]', '');
metabolite_kegg_ids = map_ids(metabolite_bigg_ids, bigg_kegg_strip_metabolite_ids(M.bigg_metabolite), M.kegg_compound);
reaction_kegg_ids   = map_ids(reaction_bigg_ids, M.bigg_reaction, M.kegg_reaction);


function out = map_ids(ids, from, to)

ind  = label_names(ids, from);
out  = repmat({''}, length(ids), 1);
out(ind ~= 0) = to(ind(ind ~= 0));
