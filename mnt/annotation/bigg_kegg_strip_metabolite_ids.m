function ids = bigg_kegg_strip_metabolite_ids(ids)

% BIGG_KEGG_STRIP_METABOLITE_IDS Remove compartment suffixes from BiGG metabolite ids (e.g. 'glc__D_e' -> 'glc')
%
% ids = bigg_kegg_strip_metabolite_ids(ids)
%
% Removes, at the end of each id (in this order)
%   1. a compartment label in brackets, e.g. '[c]', '[e]'                  ('glc[e]'  -> 'glc')
%   2. a BiGG compartment suffix, e.g. '_c', '_e', '_h', '_m'              ('13dpg_h' -> '13dpg')
%      (see bigg_metabolite_id_base for the list of compartment codes)
%   3. the stereo suffixes '__L' and '__D'                                  ('mal__L'  -> 'mal')
% as done for the ids of the IMF models. Only suffixes are removed, so '13_cis_retn_c' becomes
% '13_cis_retn'. (Until 2026-11-15, the substrings '_c', '_e', '_r', '_m', '_x', '_p', '__L', '__D'
% were removed anywhere in the id, which garbled ids such as '13_cis_retn_c' -> '13isetn'; suffixes of
% other compartments, such as '_h' (chloroplast), were kept.)
%
% See also bigg_metabolite_id_base (removes only the compartment; for matching data to models)

ids = bigg_metabolite_id_base(ids);
ids = regexprep(ids, '__[LD]$', '');
