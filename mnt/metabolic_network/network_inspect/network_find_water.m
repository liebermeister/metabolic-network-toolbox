function ind_water = network_find_water(network)

% NETWORK_FIND_WATER Find the metabolite indices of water
%
% ind_water = network_find_water(network)
%
% Water is excluded when equilibrium constants are converted between standard concentrations
% (Keq scales with 1000^(sum of stoichiometric coefficients of the dissolved compounds)), and it
% gets a fixed concentration in parameter balancing. If water is not recognised, these
% calculations are wrong by factors of 1000, so the recognition is deliberately broad.
%
% A metabolite is water if its ID or name (fields metabolites, metabolite_names, metabolite_BIGGID,
% Compound_Name), after removing an "M_" prefix and a compartment suffix ("_c", "[c]", "(c)"),
% is "water" or "h2o" (any case), or if its KEGG ID is C00001 or its ChEBI ID is 15377 (fields
% metabolite_KEGGID, Compound_Identifiers, Identifiers_kegg_compound, metabolite_ChEBI).
% Water in several compartments gives several indices.
%
% A warning is issued if a metabolite ID or name contains the word "water" or "H2O" but is not
% recognised as water (H2O2 does not count).

nm = length(network.metabolites);
is_water = false(nm,1);
looks_like_water = false(nm,1);

name_fields = {'metabolites','metabolite_names','metabolite_BIGGID','Compound_Name'};
for it = 1:length(name_fields)
  if ~isfield(network, name_fields{it}), continue; end
  names = network.(name_fields{it});
  if ~iscell(names) || length(names) ~= nm, continue; end
  for i = 1:nm
    if ~ischar(names{i}), continue; end
    n = lower(strtrim(names{i}));
    n = regexprep(n, '^m_', '');
    n = regexprep(n, '(_[a-z]{1,2}|\[[a-z0-9_]+\]|\([a-z0-9_]+\))$', '');
    if any(strcmp(n, {'water','h2o'}))
      is_water(i) = true;
    end
    if ~isempty(regexp(lower(names{i}), '(^|[^a-z0-9])(m_)?(water|h2o)([^a-z0-9]|$)', 'once'))
      looks_like_water(i) = true;
    end
  end
end

id_fields = {'metabolite_KEGGID','Compound_Identifiers','Identifiers_kegg_compound','metabolite_ChEBI'};
for it = 1:length(id_fields)
  if ~isfield(network, id_fields{it}), continue; end
  ids = network.(id_fields{it});
  if ~iscell(ids) || length(ids) ~= nm, continue; end
  for i = 1:nm
    if ~ischar(ids{i}), continue; end
    % an entry may contain several IDs, e.g. 'kegg:C00001 chebi:15377'
    tokens = regexp(lower(ids{i}), '[^\s,;|]+', 'match');
    tokens = regexprep(tokens, '^(kegg\.compound|kegg|chebi:chebi|chebi)[:_]', '');
    if any(strcmp(tokens, 'c00001')) || any(strcmp(tokens, '15377'))
      is_water(i) = true;
    end
  end
end

ind_water = find(is_water);

ind_suspicious = find(looks_like_water & ~is_water);
if length(ind_suspicious)
  warning('network_find_water: metabolite(s) %s look like water but are not recognised as water. Water must be excluded from standard concentration conversions; please check the metabolite IDs.', strjoin(column(network.metabolites(ind_suspicious))', ', '));
end
