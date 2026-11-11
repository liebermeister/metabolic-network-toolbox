function SBMLmodel = sbml_network_state_save(network, c, v, filename, dmu, options, keq)

% SBML_NETWORK_STATE_SAVE Save a network and a metabolic state (c, v, dmu, keq) to an SBML file
%
% SBMLmodel = sbml_network_state_save(network, c, v, filename, dmu, options, keq)
%
% Save a 'network' data structure, a concentration vector 'c', a flux vector 'v' and,
% optionally, reaction Gibbs free energies 'dmu' and equilibrium constants 'keq' to an
% SBML file (SBML Level 3 Version 2). Counterpart of sbml_network_state_load; analogous
% to sbtab_network_state_save.
%
% ARGUMENTS
%   network   network struct
%   c         metabolite concentrations (vector, one value per metabolite)
%   v         fluxes (vector, one value per reaction)
%   filename  name of the SBML (.xml) file; if empty, no file is written
%   dmu       (optional) reaction Gibbs free energies in kJ/mol (empty: not written)
%   options   (optional) struct with fields
%               .name               model name                          (default 'Model')
%               .concentration_unit 'M', 'mM', 'uM', 'nM'               (default 'mM')
%               .time_unit          's', 'min', 'h'                     (default 's')
%               .notes              text added to the model notes       (default '')
%               .kinetics           write the kinetic laws of network.kinetics, if present (default 1);
%                                   kinetics without field 'type' or whose sizes do not match
%                                   the network are left out (with a warning)
%               .verbose                                                (default 1)
%   keq       (optional) equilibrium constants (empty: not written)
%
% OUTPUT
%   SBMLmodel  SBML model struct in libSBML's Matlab format (can be written with OutputSBML)
%
% How the state is stored in the SBML file
%   - Concentrations: initialConcentration of the <species>, in the concentration unit
%     (substances in moles scaled to the concentration unit, e.g. mmol for mM; volumes in litres).
%   - Compartments: from network.compartments, network.compartment_sizes (in litres) and
%     network.metabolite_compartments, if present; otherwise one compartment of size 1 litre.
%   - Fluxes, equilibrium constants, reaction Gibbs free energies: local parameters of each
%     reaction (<localParameter> elements in the <kineticLaw> of the <reaction>), with the ids
%       reference_flux               flux (SBO:0000612 rate of reaction),
%                                    unit: concentration unit / time unit
%       equilibrium_constant         equilibrium constant (SBO:0000281), unit: concentration
%                                    unit ^ (sum of the stoichiometric coefficients)
%       reaction_gibbs_free_energy   Gibbs free energy of reaction (SBO:0000617), unit: kJ/mol
%     Each parameter has a unit definition and notes explaining its meaning; the model notes
%     explain the whole format.
%   - Kinetic laws: if the network has no kinetics, each kinetic law is the constant
%     rate 'reference_flux * flux_reference_volume' (amount per time), where the global
%     parameter flux_reference_volume = 1 litre is the volume to which the fluxes refer. Since the reference fluxes are stationary, a simulation of the
%     SBML model stays in the reference state. If the network has kinetics (network.kinetics),
%     the kinetic laws are written by network_sbml_export, and the three parameters are added
%     as additional local parameters (they do not appear in the rate law).
%   - External metabolites: boundaryCondition = true.
%   - Database identifiers: if the network contains the fields metabolite_KEGGID,
%     metabolite_BIGGID, metabolite_CHEBIID, reaction_KEGGID, reaction_BIGGID, the ids are
%     written as MIRIAM annotations (bqbiol:is, https://identifiers.org/kegg.compound:...,
%     bigg.metabolite, chebi, kegg.reaction, bigg.reaction; characters that are not allowed
%     in URIs are percent-encoded).
%   - Model ids: metabolite and reaction ids that are not valid SBML ids are changed by
%     network_adjust_names_for_sbml_export (with a warning); the SBML 'name' attribute
%     contains the original id.
%
% This function requires the Matlab interface of libSBML (TranslateSBML and OutputSBML)
%
% See also: sbml_network_state_load, sbtab_network_state_save, network_sbml_export

if ~exist('OutputSBML','file')
  error('Please install the Matlab interface of libSBML (https://sbml.org/software/libsbml/)');
end

eval(default('filename','[]','dmu','[]','options','struct','keq','[]'));
if isempty(options), options = struct; end
options = join_struct(struct('name','Model','concentration_unit','mM','time_unit','s','notes','','kinetics',1,'verbose',1), options);

nm = length(network.metabolites);
nr = length(network.actions);
c  = c(:);
v  = v(:);
if length(c) ~= nm, error('Length of concentration vector does not match the number of metabolites'); end
if length(v) ~= nr, error('Length of flux vector does not match the number of reactions'); end
if ~isempty(dmu), dmu = dmu(:); if length(dmu) ~= nr, error('Length of dmu vector does not match the number of reactions'); end; end
if ~isempty(keq), keq = keq(:); if length(keq) ~= nr, error('Length of keq vector does not match the number of reactions'); end; end

[conc_scale, time_multiplier] = sbml_state_units(options.concentration_unit, options.time_unit);
conc_unit_string = options.concentration_unit;
flux_unit_string = sprintf('%s/%s', options.concentration_unit, options.time_unit);

% -------------------------------------------------------------
% Network (species, reactions, kinetic laws if any) from network_sbml_export

network.s_init = c;

if isfield(network,'kinetics')
  if ~options.kinetics
    network = rmfield(network,'kinetics');
  elseif ~kinetics_sizes_ok(network.kinetics, nm, nr)
    warning('network.kinetics has no field ''type'', or its sizes do not match the network; the kinetics are not written to the SBML file');
    network = rmfield(network,'kinetics');
  else
    network.kinetics.c = c;
  end
end

[metabolites_sbml, actions_sbml] = network_adjust_names_for_sbml_export(network.metabolites, network.actions);
if ~isequal(metabolites_sbml, network.metabolites) || ~isequal(actions_sbml, network.actions)
  warning('Some metabolite or reaction ids are no valid SBML ids and were changed; the original ids are stored in the SBML name attributes');
end

SBMLmodel = network_sbml_export(network, 0, options.name, [], '', 3, 2);

if length(SBMLmodel.species) ~= nm || length(SBMLmodel.reaction) ~= nr
  error('Unexpected number of species or reactions in the exported SBML model');
end

t = sbml_state_templates();

% -------------------------------------------------------------
% Units

ud = t.unitDefinition([]);
ud(end+1) = make_unit_definition(t, 'substance_unit', {'mole', 1, conc_scale, 1});
ud(end+1) = make_unit_definition(t, 'concentration_unit', {'mole', 1, conc_scale, 1; 'litre', -1, 0, 1});
ud(end+1) = make_unit_definition(t, 'flux_unit', {'mole', 1, conc_scale, 1; 'litre', -1, 0, 1; 'second', -1, 0, time_multiplier});
if time_multiplier ~= 1
  ud(end+1) = make_unit_definition(t, 'time_unit', {'second', 1, 0, time_multiplier});
  time_unit_id = 'time_unit';
else
  time_unit_id = 'second';
end
if ~isempty(dmu)
  ud(end+1) = make_unit_definition(t, 'kJ_per_mol', {'joule', 1, 3, 1; 'mole', -1, 0, 1});
end

% units of equilibrium constants: concentration unit ^ (sum of stoichiometric coefficients)
keq_unit_id = repmat({''}, nr, 1);
n_stoich    = full(sum(network.N, 1))';
if ~isempty(keq)
  for it = 1:nr
    n = n_stoich(it);
    if n == 0
      keq_unit_id{it} = 'dimensionless';
    elseif isfinite(n)
      my_id = sprintf('concentration_unit_power_%s', strrep(strrep(num2str(n),'-','minus_'),'.','_point_'));
      if ~any(strcmp(my_id, {ud.id}))
        ud(end+1) = make_unit_definition(t, my_id, {'mole', n, conc_scale, 1; 'litre', -n, 0, 1});
      end
      keq_unit_id{it} = my_id;
    end
  end
end

SBMLmodel.unitDefinition = ud;
SBMLmodel.substanceUnits = 'substance_unit';
SBMLmodel.extentUnits    = 'substance_unit';
SBMLmodel.volumeUnits    = 'litre';
SBMLmodel.timeUnits      = time_unit_id;

% compartments

comp_template = SBMLmodel.compartment(1);
if isfield(network,'compartments') && isfield(network,'compartment_sizes') && isfield(network,'metabolite_compartments') ...
      && numel(network.metabolite_compartments) == nm && numel(network.compartment_sizes) == numel(network.compartments)
  compartment_ids = network_adjust_names_for_sbml_export(network.compartments(:));
  comp = comp_template([]);
  for it = 1:length(compartment_ids)
    comp(it)           = comp_template;
    comp(it).id        = compartment_ids{it};
    comp(it).name      = network.compartments{it};
    comp(it).size      = network.compartment_sizes(it);
    comp(it).isSetSize = double(isfinite(network.compartment_sizes(it)));
    comp(it).units     = 'litre';
  end
  SBMLmodel.compartment = comp;
  metabolite_compartment_ids = compartment_ids(label_names(network.metabolite_compartments(:), network.compartments(:)));
else
  SBMLmodel.compartment(1).units = 'litre';
  metabolite_compartment_ids = repmat({SBMLmodel.compartment(1).id}, nm, 1);
end

% global parameter: volume to which the fluxes refer

flux_volume            = t.parameter;
flux_volume.id         = 'flux_reference_volume';
flux_volume.name       = 'flux reference volume';
flux_volume.value      = 1;
flux_volume.isSetValue = 1;
flux_volume.units      = 'litre';
flux_volume.constant   = 1;
flux_volume.notes      = xhtml_notes({sprintf('Volume to which the fluxes refer (1 litre): a flux of 1 %s corresponds to a reaction rate of 1 %s per %s.', flux_unit_string, strrep(conc_unit_string,'M','mol'), options.time_unit)});
SBMLmodel.parameter = flux_volume;

% -------------------------------------------------------------
% Species: names and concentrations

metabolite_annotations = annotation_strings(network, nm, 'meta_species_', {'metabolite_KEGGID', 'kegg.compound'; 'metabolite_BIGGID', 'bigg.metabolite'; 'metabolite_CHEBIID', 'chebi'}, metabolites_sbml);
reaction_annotations   = annotation_strings(network, nr, 'meta_reaction_', {'reaction_KEGGID', 'kegg.reaction'; 'reaction_BIGGID', 'bigg.reaction'}, actions_sbml);

for it = 1:nm
  if length(metabolite_annotations{it})
    SBMLmodel.species(it).metaid     = ['meta_species_' metabolites_sbml{it}];
    SBMLmodel.species(it).annotation = metabolite_annotations{it};
  end
  SBMLmodel.species(it).name                      = network.metabolites{it};
  SBMLmodel.species(it).compartment               = metabolite_compartment_ids{it};
  SBMLmodel.species(it).initialConcentration      = c(it);
  SBMLmodel.species(it).isSetInitialConcentration = double(isfinite(c(it)));
  SBMLmodel.species(it).hasOnlySubstanceUnits     = 0;
  SBMLmodel.species(it).substanceUnits            = 'substance_unit';
end

% -------------------------------------------------------------
% Reactions: names, kinetic laws, state parameters

for it = 1:nr

  SBMLmodel.reaction(it).name = network.actions{it};
  if length(reaction_annotations{it})
    SBMLmodel.reaction(it).metaid     = ['meta_reaction_' actions_sbml{it}];
    SBMLmodel.reaction(it).annotation = reaction_annotations{it};
  end

  kl = SBMLmodel.reaction(it).kineticLaw;
  if isempty(kl)
    kl         = t.kineticLaw;
    kl.math    = 'reference_flux * flux_reference_volume';
    kl.notes   = xhtml_notes({'Constant rate: reference flux (local parameter reference_flux) times the flux reference volume (1 litre). Since the reference fluxes are stationary, a simulation keeps the model in the reference state (for compartments of size 1 litre).'});
  end

  lp = kl.localParameter;
  if isempty(lp), lp = t.localParameter([]); end

  lp(end+1) = make_local_parameter(t, 'reference_flux', v(it), 'flux_unit', 612, ...
                 sprintf('Flux (rate of reaction) in the reference state, in %s. Positive values: flux in the direction of the reaction formula (from reactants to products).', flux_unit_string));
  if ~isempty(keq)
    if isempty(keq_unit_id{it})
      unit_text = sprintf('%s to the power of the sum of the stoichiometric coefficients (products positive, reactants negative)', conc_unit_string);
    elseif strcmp(keq_unit_id{it}, 'dimensionless')
      unit_text = 'dimensionless (the sum of the stoichiometric coefficients is zero)';
    else
      unit_text = sprintf('%s^%s (sum of the stoichiometric coefficients, products positive, reactants negative)', conc_unit_string, num2str(n_stoich(it)));
    end
    lp(end+1) = make_local_parameter(t, 'equilibrium_constant', keq(it), keq_unit_id{it}, 281, ...
                 sprintf('Equilibrium constant (concentration-based, in the direction of the reaction formula: product of the equilibrium concentrations of the products divided by the product of the equilibrium concentrations of the reactants, each to the power of its stoichiometric coefficient). Unit: %s.', unit_text));
  end
  if ~isempty(dmu)
    lp(end+1) = make_local_parameter(t, 'reaction_gibbs_free_energy', dmu(it), 'kJ_per_mol', 617, ...
                 'Gibbs free energy of reaction in the reference state, in kJ/mol, in the direction of the reaction formula: Delta G = -RT ln(equilibrium_constant) + RT sum_i n_i ln(c_i). A negative value means that the reaction runs in the direction of the formula.');
  end

  kl.localParameter = lp;
  SBMLmodel.reaction(it).kineticLaw = kl;

end

% -------------------------------------------------------------
% Model notes

notes_text = {sprintf('Model %s with a metabolic reference state, written by sbml_network_state_save (Metabolic Network Toolbox, https://github.com/liebermeister/mnt).', options.name), ...
              sprintf('Metabolite concentrations: initial concentrations of the species, in %s (compartment sizes in litres).', conc_unit_string), ...
              sprintf('Fluxes: local parameter "reference_flux" of each reaction (SBO:0000612 rate of reaction), in %s.', flux_unit_string)};
if ~isempty(keq)
  notes_text{end+1} = sprintf('Equilibrium constants: local parameter "equilibrium_constant" of each reaction (SBO:0000281), concentration-based with concentrations in %s.', conc_unit_string);
end
if ~isempty(dmu)
  notes_text{end+1} = 'Reaction Gibbs free energies: local parameter "reaction_gibbs_free_energy" of each reaction (SBO:0000617), in kJ/mol.';
end
notes_text{end+1} = 'All reaction quantities refer to the direction of the reaction formula (reactants to products). External metabolites have boundaryCondition = true.';
if isfield(network,'kinetics')
  notes_text{end+1} = 'The kinetic laws describe the model kinetics; the local parameters reference_flux, equilibrium_constant and reaction_gibbs_free_energy are not used in them.';
else
  notes_text{end+1} = 'Kinetic laws: constant rates, reference flux times flux_reference_volume (1 litre); the model contains no other kinetic information.';
end
if length(options.notes)
  notes_text{end+1} = options.notes;
end
SBMLmodel.notes = xhtml_notes(notes_text);

% -------------------------------------------------------------
% Write file

if length(filename)
  OutputSBML(SBMLmodel, filename);
  if options.verbose
    fprintf('Wrote SBML file "%s"\n', filename);
  end
end


% -------------------------------------------------------------------

function annotations = annotation_strings(network, n, metaid_prefix, field_collection, ids)

% MIRIAM annotations (RDF) with database identifiers from network fields
% field_collection: one row per database: {network field, identifiers.org collection}

annotations = repmat({''}, n, 1);
for it = 1:n
  resources = {};
  for itf = 1:size(field_collection,1)
    fn = field_collection{itf,1};
    if isfield(network, fn) && numel(network.(fn)) == n
      my_id = network.(fn){it};
      if ischar(my_id) && length(strtrim(my_id))
        resources{end+1} = sprintf('<rdf:li rdf:resource="https://identifiers.org/%s:%s"/>', field_collection{itf,2}, uri_encode(strtrim(my_id)));
      end
    end
  end
  if length(resources)
    annotations{it} = ['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">' ...
                       '<rdf:Description rdf:about="#' metaid_prefix ids{it} '"><bqbiol:is><rdf:Bag>' strjoin(resources,'') ...
                       '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>'];
  end
end


% -------------------------------------------------------------------

function s = uri_encode(s)

% percent-encode all characters except letters, digits and - . _ ~ :

ind = find(~(isletter(s) & s < 128) & ~ismember(s, '0123456789-._~:'));
for it = fliplr(ind)
  s = [s(1:it-1) sprintf('%%%02X', double(s(it))) s(it+1:end)];
end


% -------------------------------------------------------------------

function ok = kinetics_sizes_ok(kinetics, nm, nr)

% check the sizes of the kinetic parameters used by network_sbml_export

ok = 1;
if ~isfield(kinetics,'type'), ok = 0; return; end
switch kinetics.type
  case {'ms','cs','ds','rp'}
    for fn = {'u','KV','Keq'}
      if isfield(kinetics,fn{1}) && ~isempty(kinetics.(fn{1})) && numel(kinetics.(fn{1})) ~= nr, ok = 0; end
    end
    if ~isfield(kinetics,'KV') || ~isfield(kinetics,'Keq'), ok = 0; end
    for fn = {'KM','KA','KI'}
      if ~isfield(kinetics,fn{1}) || ~isequal(size(kinetics.(fn{1})), [nr, nm]), ok = 0; end
    end
  case 'kinetic_strings'
    if ~isfield(kinetics,'reactions') || numel(kinetics.reactions) ~= nr, ok = 0; end
end


% -------------------------------------------------------------------

function [conc_scale, time_multiplier] = sbml_state_units(concentration_unit, time_unit)

switch concentration_unit
  case 'M',  conc_scale =  0;
  case 'mM', conc_scale = -3;
  case 'uM', conc_scale = -6;
  case 'nM', conc_scale = -9;
  otherwise, error('Unknown concentration unit "%s" (allowed: M, mM, uM, nM)', concentration_unit);
end

switch time_unit
  case 's',   time_multiplier = 1;
  case 'min', time_multiplier = 60;
  case 'h',   time_multiplier = 3600;
  otherwise, error('Unknown time unit "%s" (allowed: s, min, h)', time_unit);
end


% -------------------------------------------------------------------

function t = sbml_state_templates()

% Empty libSBML structs (SBML Level 3 Version 2) for unit definitions, units, kinetic laws
% and local parameters, obtained by reading a minimal SBML file with TranslateSBML

xml = ['<?xml version="1.0" encoding="UTF-8"?>' newline ...
       '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">' newline ...
       '<model id="m">' newline ...
       '<listOfUnitDefinitions><unitDefinition id="u"><listOfUnits><unit kind="mole" exponent="1" scale="0" multiplier="1"/></listOfUnits></unitDefinition></listOfUnitDefinitions>' newline ...
       '<listOfCompartments><compartment id="c" constant="true"/></listOfCompartments>' newline ...
       '<listOfSpecies><species id="A" compartment="c" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>' newline ...
       '<listOfParameters><parameter id="p" constant="true"/></listOfParameters>' newline ...
       '<listOfReactions><reaction id="r" reversible="true">' newline ...
       '<listOfReactants><speciesReference species="A" constant="true"/></listOfReactants>' newline ...
       '<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>k</ci></math>' newline ...
       '<listOfLocalParameters><localParameter id="k"/></listOfLocalParameters></kineticLaw>' newline ...
       '</reaction></listOfReactions>' newline ...
       '</model>' newline '</sbml>' newline];

filename = [tempname '.xml'];
fid = fopen(filename,'w'); fprintf(fid,'%s',xml); fclose(fid);
m = TranslateSBML(filename, 0, 0);
delete(filename);

t.parameter                 = m.parameter(1);
t.unit                      = m.unitDefinition(1).unit;
t.unitDefinition            = m.unitDefinition(1);
t.unitDefinition.unit       = t.unit([]);
t.localParameter            = m.reaction(1).kineticLaw.localParameter;
t.kineticLaw                = m.reaction(1).kineticLaw;
t.kineticLaw.localParameter = t.localParameter([]);


% -------------------------------------------------------------------

function ud = make_unit_definition(t, id, units)

% units: cell array, one row per unit: {kind, exponent, scale, multiplier}

ud    = t.unitDefinition;
ud.id = id;
for it = 1:size(units,1)
  u            = t.unit;
  u.kind       = units{it,1};
  u.exponent   = units{it,2};
  u.scale      = units{it,3};
  u.multiplier = units{it,4};
  ud.unit(it)  = u;
end


% -------------------------------------------------------------------

function p = make_local_parameter(t, id, value, units, sbo_term, note)

p            = t.localParameter;
p.id         = id;
p.name       = strrep(id,'_',' ');
p.value      = value;
p.isSetValue = double(isfinite(value));
p.units      = units;
p.sboTerm    = sbo_term;
p.notes      = xhtml_notes({note});


% -------------------------------------------------------------------

function notes = xhtml_notes(paragraphs)

notes = '<notes><body xmlns="http://www.w3.org/1999/xhtml">';
for it = 1:length(paragraphs)
  notes = [notes '<p>' xml_escape(paragraphs{it}) '</p>'];
end
notes = [notes '</body></notes>'];


% -------------------------------------------------------------------

function s = xml_escape(s)

s = strrep(s,'&','&amp;');
s = strrep(s,'<','&lt;');
s = strrep(s,'>','&gt;');
