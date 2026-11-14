% Tests for sbml_network_state_save and sbml_network_state_load (requires libSBML's Matlab interface)

setup_test_path;
if ~exist('TranslateSBML', 'file')
  disp('libSBML not installed; skipping SBML state test');
  return
end

network = example_network;
network.metabolites = {'S1', 'S2', 'S3[c]', 'S4', 'S5', 'S6'}';   % 'S3[c]' is no valid SBML id
network.metabolite_KEGGID = {'C00001', '', 'C00003', '', '', ''}';
network.reaction_BIGGID   = {'R1', 'R2', 'R3', 'R4', 'R5'}';
network.compartments            = {'cell', 'medium'};
network.compartment_sizes       = [1 2];
network.metabolite_compartments = {'medium', 'medium', 'cell', 'cell', 'medium', 'medium'}';
c   = [2 1 0.5 0.4 0.1 0.2]';
v   = [1 1 2 1.5 0.5]';                       % stationary: N(internal,:) * v = 0
keq = [10 20 5 3 4]';
dmu = -RT * (log(keq) - network.N' * log(c));

%% file round trip
filename = [tempname '.xml'];
warning('off'); 
sbml_network_state_save(network, c, v, filename, dmu, struct('verbose', 0), keq);
[network2, c2, v2, dmu2, keq2, theta2] = sbml_network_state_load(filename);
warning('on');
[~, err] = TranslateSBML(filename, 1, 0);
delete(filename);
assert(isempty(err));                                  % valid SBML, no warnings
assert(isequal(network2.metabolites(:), network.metabolites(:)));   % original ids restored
assert(isequal(network2.actions(:), network.actions(:)));
assert(isequal(full(network2.N), full(network.N)));
assert(isequal(network2.external(:), network.external(:)));
assert(isequal(network2.metabolite_compartments(:), network.metabolite_compartments(:)));
assert(isequal(network2.compartment_sizes(:), network.compartment_sizes(:)));
assert(isequal(network2.metabolite_KEGGID(:), network.metabolite_KEGGID(:)));
assert(isequal(network2.reaction_BIGGID(:), network.reaction_BIGGID(:)));
assert(~isfield(network2, 'kinetics'));
assert(max(abs(c2 - c)) < 1e-12 && max(abs(v2 - v)) < 1e-12);
assert(max(abs(keq2 - keq)) < 1e-12 && max(abs(dmu2 - dmu)) < 1e-10);
assert(max(abs(theta2 - (log(keq) - network.N' * log(c)))) < 1e-12);

%% in memory, without keq and dmu, other units
sbml = sbml_network_state_save(network, c, v, [], [], struct('concentration_unit', 'uM', 'time_unit', 'min', 'verbose', 0));
warning('off'); [~, c3, v3, dmu3, keq3] = sbml_network_state_load(sbml); warning('on');
assert(isempty(dmu3) && isempty(keq3));
assert(max(abs(c3 - 0.001 * c)) < 1e-15 && max(abs(v3 - v)) < 1e-12);   % values written as uM, returned in mM

%% file in M: concentrations returned in mM, equilibrium constants converted to 1 mM
sbml = sbml_network_state_save(network, c, v, [], [], struct('concentration_unit', 'M', 'verbose', 0), keq);
[~, c5, ~, ~, keq5] = sbml_network_state_load(sbml);
Nw = network.N; Nw(network_find_water(network),:) = 0;
assert(max(abs(c5 - 1000 * c)) < 1e-9);
assert(max(abs(keq5 ./ (keq .* 1000 .^ full(sum(Nw,1))') - 1)) < 1e-12);

%% network with kinetics: rate laws kept, state parameters removed from the kinetics
network.kinetics = set_kinetics(network, 'cs');
sbml = sbml_network_state_save(network, c, v, [], dmu, struct('verbose', 0), keq);
warning('off'); [network4, c4, v4, ~, keq4] = sbml_network_state_load(sbml); warning('on');
assert(strcmp(network4.kinetics.type, 'kinetic_strings'));
assert(~any(ismember(network4.kinetics.reactions{1}.parameters, {'reference_flux', 'equilibrium_constant', 'reaction_gibbs_free_energy'})));
assert(max(abs(v4 - v)) < 1e-12 && max(abs(keq4 - keq)) < 1e-12);
