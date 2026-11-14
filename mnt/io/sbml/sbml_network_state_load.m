function [network, c, v, dmu, keq, theta] = sbml_network_state_load(filename)

% SBML_NETWORK_STATE_LOAD Load a network and a metabolic state (c, v, dmu, keq) from an SBML file
%
% [network, c, v, dmu, keq, theta] = sbml_network_state_load(filename)
%
% Load a 'network' data structure, a concentration vector 'c', and a flux vector 'v' from an
% SBML file written by sbml_network_state_save (and, optionally, reaction Gibbs free energies
% 'dmu' and equilibrium constants 'keq'; empty if not in the file). Analogous to
% sbtab_network_state_load.
%
% filename: name of an SBML file, or an SBML model struct (as returned by TranslateSBML)
%
% Where the values are read from (see sbml_network_state_save)
%   c    initial concentrations of the species
%   v    local parameter 'reference_flux' of each reaction
%   keq  local parameter 'equilibrium_constant' of each reaction
%   dmu  local parameter 'reaction_gibbs_free_energy' of each reaction
%   theta = log(keq) - network.N' * log(c)   (thermodynamic forces -dmu/RT; empty if keq is empty)
% For files written by sbml_network_state_save, the metabolite and reaction ids are taken
% from the SBML name attributes (original ids, also if they are no valid SBML ids).
% Database identifiers in MIRIAM annotations (identifiers.org: kegg.compound, bigg.metabolite,
% chebi, kegg.reaction, bigg.reaction) are stored in the network fields metabolite_KEGGID,
% metabolite_BIGGID, metabolite_CHEBIID, reaction_KEGGID, reaction_BIGGID.
% Missing values are NaN. Concentrations and equilibrium constants are converted to mM (Keq: standard
% concentration 1 mM) from the substance unit declared in the file (unit definition of the model's
% substanceUnits, kind 'mole'); without such a declaration, mM is assumed (with a warning). Fluxes and
% Gibbs free energies are returned as stored (default: mM/s, kJ/mol).
%
% If the kinetic laws are only the constant rates 'reference_flux * flux_reference_volume' (files written for networks
% without kinetics), the network has no field 'kinetics'. Otherwise, the kinetics are read by
% network_sbml_import, without the three state parameters.
%
% This function requires the Matlab interface of libSBML (TranslateSBML)
%
% See also: sbml_network_state_save, sbtab_network_state_load, network_sbml_import

if ~exist('TranslateSBML','file')
  error('Please install the Matlab interface of libSBML (https://sbml.org/software/libsbml/)');
end

if ischar(filename) || isstring(filename)
  if ~exist(filename,'file'), error('File "%s" not found', filename); end
  s = TranslateSBML(char(filename), 0, 0);
else
  s = filename;
end

network = network_sbml_import(s);

nm = length(s.species);
nr = length(s.reaction);

% -------------------------------------------------------------
% Original metabolite and reaction ids (files written by sbml_network_state_save store them
% in the name attributes, because ids that are no valid SBML ids are changed on export)

written_by_state_save = isfield(s,'parameter') && any(strcmp('flux_reference_volume', {s.parameter.id}));
if written_by_state_save
  species_names  = {s.species.name}';
  reaction_names = {s.reaction.name}';
  if all(~cellfun(@isempty, species_names)) && length(unique(species_names)) == nm
    network.metabolites      = species_names;
  end
  if all(~cellfun(@isempty, reaction_names)) && length(unique(reaction_names)) == nr
    network.actions          = reaction_names;
  end
end

% -------------------------------------------------------------
% Database identifiers from MIRIAM annotations (identifiers.org URIs)

network = read_identifiers(network, {s.species.annotation},  {'kegg.compound', 'metabolite_KEGGID'; 'bigg.metabolite', 'metabolite_BIGGID'; 'chebi', 'metabolite_CHEBIID'});
network = read_identifiers(network, {s.reaction.annotation}, {'kegg.reaction', 'reaction_KEGGID'; 'bigg.reaction', 'reaction_BIGGID'});

% -------------------------------------------------------------
% Concentrations

c = nan(nm,1);
for it = 1:nm
  if s.species(it).isSetInitialConcentration
    c(it) = s.species(it).initialConcentration;
  elseif isfield(s.species(it),'isSetInitialAmount') && s.species(it).isSetInitialAmount
    c(it) = s.species(it).initialAmount;
  end
end

% -------------------------------------------------------------
% Reaction parameters

state_parameter_ids = {'reference_flux', 'equilibrium_constant', 'reaction_gibbs_free_energy'};

v   = nan(nr,1);
keq = nan(nr,1);
dmu = nan(nr,1);
only_flux_laws = 1;

for it = 1:nr
  kl = s.reaction(it).kineticLaw;
  if isempty(kl)
    only_flux_laws = 0;
    continue
  end
  if isfield(kl,'localParameter'), lp = kl.localParameter; else, lp = kl.parameter; end
  for itt = 1:length(lp)
    if lp(itt).isSetValue
      switch lp(itt).id
        case 'reference_flux',             v(it)   = lp(itt).value;
        case 'equilibrium_constant',       keq(it) = lp(itt).value;
        case 'reaction_gibbs_free_energy', dmu(it) = lp(itt).value;
      end
    end
  end
  if ~any(strcmp(regexprep(kl.math,'\s',''), {'reference_flux', 'reference_flux*flux_reference_volume'}))
    only_flux_laws = 0;
  end
end

if all(isnan(v))
  warning('No flux values (local parameters reference_flux) found in file');
end

if all(isnan(dmu))
  dmu = [];
  if nargout > 3, warning('No reaction Gibbs free energy values found in file'); end
end

% -------------------------------------------------------------
% Units: concentrations and equilibrium constants are converted to mM (Keq: to the standard
% concentration 1 mM) if the model declares its substance unit (unit definition with kind 'mole')
% and the volume unit litre; otherwise the values are returned as stored (mM assumed)

[conc_factor, unit_found] = sbml_concentration_factor(s);
if unit_found && conc_factor ~= 1
  c   = conc_factor * c;
  N = network.N; N(network_find_water(network),:) = 0;   % water: activity 1, not converted
  m   = full(sum(N, 1))';
  keq = keq .* conc_factor .^ m;
  display(sprintf('- (sbml_network_state_load.m): converting concentrations and equilibrium constants from the substance unit of the file to mM (factor %g)', conc_factor));
elseif ~unit_found && (any(isfinite(c)) || any(isfinite(keq)))
  warning('sbml_network_state_load:units', 'SBML file without substance unit definition: concentrations and equilibrium constants are assumed to be in mM (standard concentration 1 mM)');
end

if all(isnan(keq))
  keq = [];
  if nargout > 4, warning('No equilibrium constant values found in file'); end
end

% -------------------------------------------------------------
% Kinetics: remove constant-flux laws, or remove the state parameters from the kinetics

if isfield(network,'kinetics')
  if only_flux_laws
    network = rmfield(network,'kinetics');
  elseif strcmp(network.kinetics.type,'kinetic_strings')
    for it = 1:length(network.kinetics.reactions)
      r = network.kinetics.reactions{it};
      ind_keep = ~ismember(r.parameters, state_parameter_ids);
      r.parameters       = r.parameters(ind_keep);
      r.parameter_values = r.parameter_values(ind_keep);
      network.kinetics.reactions{it} = r;
    end
  end
end

theta = [];
if nargout > 5 && ~isempty(keq)
  theta = log(keq) - network.N' * log(c);
end


% -------------------------------------------------------------------

function network = read_identifiers(network, annotations, collection_field)

% collection_field: one row per database: {identifiers.org collection, network field}

for itc = 1:size(collection_field,1)
  pattern = ['identifiers\.org/' strrep(collection_field{itc,1},'.','\.') '[:/]([^"<\s]+)'];
  ids   = repmat({''}, length(annotations), 1);
  found = 0;
  for it = 1:length(annotations)
    tok = regexp(annotations{it}, pattern, 'tokens', 'once');
    if ~isempty(tok)
      ids{it} = uri_decode(tok{1});
      found   = 1;
    end
  end
  if found
    network.(collection_field{itc,2}) = ids;
  end
end


% -------------------------------------------------------------------

function s = uri_decode(s)

[tok, split] = regexp(s, '%([0-9A-Fa-f]{2})', 'tokens', 'split');
out = split{1};
for it = 1:length(tok)
  out = [out char(hex2dec(tok{it}{1})) split{it+1}];
end
s = out;


% -------------------------------------------------------------

function [f, found] = sbml_concentration_factor(s)
% factor from the model's substance unit (per litre) to mM; found = 0 if no such unit definition
f = 1; found = 0;
if ~isfield(s,'substanceUnits') || isempty(s.substanceUnits) || ~isfield(s,'unitDefinition'), return; end
ud = s.unitDefinition;
k = find(strcmp({ud.id}, s.substanceUnits));
if isempty(k), return; end
u = ud(k(1)).unit;
im = find(strcmp({u.kind}, 'mole'));
if length(u) ~= 1 || length(im) ~= 1 || u(im).exponent ~= 1, return; end
mult = 1; if isfield(u(im),'multiplier') && ~isempty(u(im).multiplier), mult = double(u(im).multiplier); end
f = mult * 10^(double(u(im).scale) + 3);   % TranslateSBML returns integer scales
found = 1;
