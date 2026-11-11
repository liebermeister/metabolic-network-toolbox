function ids = bigg_kegg_strip_metabolite_ids(ids)

% BIGG_KEGG_STRIP_METABOLITE_IDS Remove compartment suffixes from BiGG metabolite ids (e.g. 'glc__D_e' -> 'glc')
%
% ids = bigg_kegg_strip_metabolite_ids(ids)
%
% Removes the substrings '_c', '_e', '_r', '_m', '_x', '_p', '__L', '__D' (in this order), as done for
% the ids of the IMF models.

for s = {'_c', '_e', '_r', '_m', '_x', '_p', '__L', '__D'}
  ids = strrep(ids, s{1}, '');
end
