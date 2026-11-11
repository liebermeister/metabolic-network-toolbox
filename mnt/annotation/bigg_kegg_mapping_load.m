function M = bigg_kegg_mapping_load(id_mapping_file)

% BIGG_KEGG_MAPPING_LOAD Load a BiGG-KEGG mapping table (e.g. bigg-kegg-mapping.tsv in this folder)
%
% M = bigg_kegg_mapping_load(id_mapping_file)
%
% The table has a header row with the columns !Identifiers:bigg.reaction, !Identifiers:bigg.metabolite,
% !Identifiers:kegg.reaction, !Identifiers:kegg.compound (other columns are ignored).
% Output: struct with the column cell arrays bigg_reaction, bigg_metabolite, kegg_reaction, kegg_compound

if ~exist('id_mapping_file', 'var') || isempty(id_mapping_file)
  id_mapping_file = fullfile(fileparts(mfilename('fullpath')), 'bigg-kegg-mapping.tsv');
end

T = load_any_table(id_mapping_file);
header_row = find(strncmp(T(:,1), '!', 1) & ~strncmp(T(:,1), '!!', 2), 1);
if isempty(header_row), error('No header row found in %s', id_mapping_file); end
header = T(header_row, :);
T = T(header_row+1:end, :);

names  = {'bigg_reaction', 'bigg_metabolite', 'kegg_reaction', 'kegg_compound'};
titles = {'!Identifiers:bigg.reaction', '!Identifiers:bigg.metabolite', '!Identifiers:kegg.reaction', '!Identifiers:kegg.compound'};
for it = 1:length(names)
  ind = find(strcmp(header, titles{it}), 1);
  if isempty(ind), error('Column %s not found in %s', titles{it}, id_mapping_file); end
  M.(names{it}) = T(:, ind);
end
