% Tests for sbtab_to_sbml and sbml_to_sbtab with metabolic states (requires libSBML's Matlab interface)

setup_test_path;
if ~exist('TranslateSBML', 'file')
  disp('libSBML not installed; skipping SBML-SBtab conversion test');
  return
end

network = example_network;
network.metabolite_KEGGID = {'C00001', '', '', '', '', ''}';   % S1: water
c   = [2 1 0.5 0.4 0.1 0.2]';                                   % mM
v   = [1 1 2 1.5 0.5]';                                         % mM/s
keq = [10 20 5 3 4]';                                           % standard concentration 1 mM
dmu = -RT * (log(keq) - network.N' * log(c));
Nw  = network.N; Nw(network_find_water(network),:) = 0;
m   = full(sum(Nw, 1))';

sbtab = network_to_sbtab(network, struct('c', c, 'v', v, 'dmu', dmu, 'keq', keq, 'verbose', 0));
% (variables defined here are available in all test sections)
warning('off'); sbml = sbtab_to_sbml(sbtab, [], struct('concentration_unit', 'uM', 'time_unit', 'min')); warning('on');

%% SBtab -> SBML in uM and min: values converted to the units of the SBML file
assert(max(abs([sbml.species.initialConcentration]' - 1000 * c)) < 1e-9);
v_sbml   = nan(5,1);
keq_sbml = nan(5,1);
for it = 1:5
  lp = sbml.reaction(it).kineticLaw.localParameter;
  v_sbml(it)   = lp(strcmp({lp.id}, 'reference_flux')).value;
  keq_sbml(it) = lp(strcmp({lp.id}, 'equilibrium_constant')).value;
end
assert(max(abs(v_sbml - 60000 * v)) < 1e-6);                        % mM/s -> uM/min
assert(max(abs(keq_sbml ./ (keq .* 1000 .^ m) - 1)) < 1e-12);       % 1 mM -> uM (water not counted)

%% SBML -> SBtab: values back in mM, mM/s, Keq for 1 mM
warning('off');
sbtab2 = sbml_to_sbtab(sbml);
[~, c2, v2, dmu2, keq2] = sbtab_network_state_load(sbtab2);
warning('on');
assert(max(abs(c2 - c)) < 1e-12);
assert(max(abs(v2 - v)) < 1e-12);
assert(max(abs(keq2 ./ keq - 1)) < 1e-12);
assert(max(abs(dmu2 - dmu)) < 1e-10);

%% option metabolite_state = 0: only the network
warning('off'); sbtab3 = sbml_to_sbtab(sbml, [], struct('metabolite_state', 0)); warning('on');
assert(~any(ismember({'Flux', 'MetaboliteConcentration'}, sbtab_document_get_table_names(sbtab3))));
warning('off'); sbml4 = sbtab_to_sbml(sbtab, [], struct('metabolite_state', 0)); warning('on');
assert(~any(strcmp('flux_reference_volume', {sbml4.parameter.id})));

%% older option name 'state'
warning('off'); sbml5 = sbtab_to_sbml(sbtab, [], struct('state', 0)); warning('on');
assert(~any(strcmp('flux_reference_volume', {sbml5.parameter.id})));

%% flux unit declared in the SBtab table Flux (table attribute Unit)
sbtab6 = network_to_sbtab(network, struct('c', c, 'v', 60000 * v, 'verbose', 0));
sbtab6.tables.Flux = sbtab_table_add_attribute(sbtab6.tables.Flux, 'Unit', 'uM/min');
warning('off'); [~, ~, v6] = sbtab_network_state_load(sbtab6); warning('on');
assert(max(abs(v6 - v)) < 1e-12);
assert(abs(flux_unit_in_mM_per_s('uM/min') - 1/60000) < 1e-20);

%% tables are found by their TableType (e.g. TableIDs 'reaction', 'compound' in files written by Python sbml2sbtab)
sbtab7 = sbtab;
sbtab7.tables.reaction = sbtab7.tables.Reaction;  sbtab7.tables = rmfield(sbtab7.tables, 'Reaction');
sbtab7.tables.compound = sbtab7.tables.Compound;  sbtab7.tables = rmfield(sbtab7.tables, 'Compound');
sbtab7.table_names = fieldnames(sbtab7.tables);
network7 = sbtab_to_network(sbtab7, struct('verbose', 0, 'omit_water', 0));
assert(isequal(network7.actions(:), network.actions(:)));
assert(isequal(network7.metabolites(:), network.metabolites(:)));
assert(isequal(full(network7.N), full(network.N)));
