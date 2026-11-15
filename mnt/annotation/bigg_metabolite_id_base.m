function [base, compartment] = bigg_metabolite_id_base(ids)

% BIGG_METABOLITE_ID_BASE Split BiGG-style metabolite ids into base id and compartment
%
% [base, compartment] = bigg_metabolite_id_base(ids)
%
% ids          metabolite id (string) or cell array of ids
% base         ids without compartment (same type as ids)
% compartment  compartment code ('' if the id has none)
%
% Two notations of compartments are recognised (only at the end of an id):
%   - label in brackets:  'glc[e]'  -> base 'glc',   compartment 'e'
%   - BiGG suffix:        'glc__D_e' -> base 'glc__D', compartment 'e';  '13dpg_h' -> '13dpg', 'h'
%     The suffix must consist of one underscore and a BiGG compartment code
%       c e p m x r v n g u l h f s i w y im cx um cm mm
%     and must not follow another underscore (stereo suffixes such as 'mal__L' are kept: base 'mal__L').
% Stereo suffixes ('__L', '__D') are part of the base id (to remove them, see bigg_kegg_strip_metabolite_ids).
%
% Used by bigg_match_metabolite_ids to compare metabolite ids of data sets and models.

compartment_codes = 'c|e|p|m|x|r|v|n|g|u|l|h|f|s|i|w|y|im|cx|um|cm|mm';

was_char = ischar(ids);
if was_char, ids = {ids}; end

base        = ids;
compartment = repmat({''}, size(ids));

for it = 1:numel(ids)
  id = ids{it};
  if ~ischar(id) || isempty(id), continue; end
  t = regexp(id, '^(.+)\[([A-Za-z0-9_]+)\]$', 'tokens', 'once');
  if isempty(t)
    t = regexp(id, ['^(.*[^_])_(' compartment_codes ')$'], 'tokens', 'once');
  end
  if ~isempty(t)
    base{it}        = t{1};
    compartment{it} = t{2};
  end
end

if was_char
  base        = base{1};
  compartment = compartment{1};
end
