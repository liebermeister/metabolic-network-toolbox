% Tests for network_find_water: water is recognised under the usual IDs and names;
% hydrogen peroxide is not water

setup_test_path;

% IDs that must be recognised as water
ids_water = {'H2O', 'h2o', 'water', 'Water', 'h2o_c', 'M_h2o_c', 'h2o_e', 'H2O[c]', 'h2o(e)'};
for it = 1:length(ids_water)
  network = struct('metabolites', {{'A'; ids_water{it}; 'B'}});
  [~, ind] = evalc('network_find_water(network)');
  assert(isequal(ind, 2), sprintf('water not recognised: %s', ids_water{it}));
end

% IDs that are not water
ids_other = {'H2O2', 'h2o2_c', 'A', 'hydroxide'};
for it = 1:length(ids_other)
  network = struct('metabolites', {{'A'; ids_other{it}; 'B'}});
  [~, ind] = evalc('network_find_water(network)');
  assert(~any(ind == 2), sprintf('wrongly recognised as water: %s', ids_other{it}));
end

% recognition by KEGG or ChEBI ID, with an arbitrary metabolite ID
network = struct('metabolites', {{'A'; 'x1'; 'x2'}}, 'metabolite_KEGGID', {{'C00002'; 'C00001'; ''}});
assert(isequal(network_find_water(network), 2));
network = struct('metabolites', {{'A'; 'x1'; 'x2'}}, 'Compound_Identifiers', {{''; ''; 'kegg:C00001 chebi:15377'}});
assert(isequal(network_find_water(network), 3));

% water in two compartments
network = struct('metabolites', {{'h2o_c'; 'A'; 'h2o_e'}});
assert(isequal(network_find_water(network), [1; 3]));

% a metabolite that looks like water but is not recognised gives a warning
network = struct('metabolites', {{'A'; 'H2O_cytosol_pool'}});
lastwarn('');
evalc('network_find_water(network)');
[msg, ~] = lastwarn;
assert(~isempty(strfind(msg, 'look like water')));
