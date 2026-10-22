% Tests for SBML export and import (requires libSBML's Matlab interface)
%
% The model is converted in memory and not written to a file

setup_test_path;
network = example_network;
network.kinetics = set_kinetics(network, 'cs');

%% network_sbml_export and network_sbml_import: round trip
if ~exist('TranslateSBML', 'file')
  disp('libSBML not installed; skipping SBML test');
  return
end
for level_version = [2 2; 2 4; 3 1]'
  sbml_model = network_sbml_export(network, 0, 'test model', [], '', level_version(1), level_version(2));
  assert(numel(sbml_model.species) == 6 && numel(sbml_model.reaction) == 5);
  network2 = network_sbml_import(sbml_model);
  assert(isequal(full(network2.N), network.N));
  assert(isequal(network2.reversible(:), network.reversible(:)));
end
