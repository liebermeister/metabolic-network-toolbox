function sbtab_object = sbml_to_sbtab(sbml_object, filename, options)

% sbtab_object = sbml_to_sbtab(sbml_object, filename, options)
%
% Convert an SBML model (libSBML Matlab struct, or SBML file name) into an SBtab document.
% If 'filename' is given, the SBtab document is also written to a file.
%
% Metabolic state: if the SBML model contains the local reaction parameters 'reference_flux',
% 'equilibrium_constant' or 'reaction_gibbs_free_energy' (SBML format of the Metabolic Network
% Toolbox for models with a metabolic state, written by sbml_network_state_save or by the Python
% converter sbtab2sbml --metabolite-state), the model is read with sbml_network_state_load, and the
% SBtab document contains the tables MetaboliteConcentration (from the initial concentrations of
% the species, in mM), Flux (in mM/s), ReactionGibbsFreeEnergy (in kJ/mol) and EquilibriumConstant
% (standard concentration 1 mM), as written by sbtab_network_state_save. Otherwise, only the
% network is converted (network_sbml_import, network_to_sbtab).
%
% options (optional struct)
%   .metabolite_state   'auto' (default: read the state if the SBML model contains it),
%                       1 (always read the state format), 0 (never: only the network is converted)
%
% This function requires the Metabolic Network Toolbox and libSBML (TranslateSBML)
%
% usage example:
% sbml_object = sbtab_load_example('sbml/escherichia-coli_noor_2016.xml','sbml');
% sbtab_object = sbml_to_sbtab(sbml_object);
%
% See also: sbtab_to_sbml, sbml_network_state_load, sbtab_network_state_save

eval(default('filename','[]','options','struct'));
if isempty(options), options = struct; end
options = join_struct(struct('metabolite_state','auto'), options);

if ischar(sbml_object) || isstring(sbml_object)
  sbml_object = TranslateSBML(char(sbml_object), 0, 0);
end

if ischar(options.metabolite_state) && strcmp(options.metabolite_state, 'auto')
  read_state = sbml_contains_state(sbml_object);
else
  read_state = options.metabolite_state;
end

if read_state
  warning_state = warning('off', 'all');   % missing state values are fine here
  [network, c, v, dmu, keq] = sbml_network_state_load(sbml_object);
  warning(warning_state);
  % fluxes: from the unit of the local parameters reference_flux to mM/s
  f_v = sbml_flux_factor(sbml_object);
  if isfinite(f_v)
    v = f_v * v;
  elseif any(isfinite(v))
    warning('Unknown unit of the fluxes (local parameters reference_flux) in the SBML model; mM/s is assumed');
  end
  if all(isnan(c)), c = []; end
  if all(isnan(v)), v = []; end
  sbtab_object = network_to_sbtab(network, struct('c', c, 'v', v, 'dmu', dmu, 'keq', keq, 'verbose', 0));
else
  sbtab_object = network_to_sbtab(network_sbml_import(sbml_object), struct('verbose', 0));
end

if length(filename)
  sbtab_document_save(sbtab_object, filename, 0, 0);
end


% -------------------------------------------------------------------

function flag = sbml_contains_state(s)

flag = 0;
for it = 1:length(s.reaction)
  lp = reaction_local_parameters(s.reaction(it));
  if ~isempty(lp) && any(ismember({lp.id}, {'reference_flux', 'equilibrium_constant', 'reaction_gibbs_free_energy'}))
    flag = 1;
    return
  end
end


% -------------------------------------------------------------------

function lp = reaction_local_parameters(reaction)

lp = [];
kl = reaction.kineticLaw;
if isempty(kl), return; end
if isfield(kl,'localParameter'), lp = kl.localParameter; else, lp = kl.parameter; end


% -------------------------------------------------------------------

function f = sbml_flux_factor(s)

% Factor from the unit of the local parameters reference_flux to mM/s (NaN if unknown).
% Without unit (e.g. files without unit definitions), mM/s is assumed (factor 1).

units = {};
for it = 1:length(s.reaction)
  lp = reaction_local_parameters(s.reaction(it));
  for itt = 1:length(lp)
    if strcmp(lp(itt).id, 'reference_flux') && isfield(lp(itt), 'units')
      units{end+1} = lp(itt).units;
    end
  end
end
units = unique(units);
f = 1;
if isempty(units) || (length(units) == 1 && isempty(units{1})), return; end
if length(units) > 1, f = nan; return; end
f = nan;
if ~isfield(s, 'unitDefinition') || isempty(s.unitDefinition), return; end
k = find(strcmp({s.unitDefinition.id}, units{1}));
if isempty(k), return; end
u = s.unitDefinition(k(1)).unit;
% the unit must be (amount of substance) / (volume * time)
expected = struct('mole', 1, 'litre', -1, 'second', -1);
if length(u) ~= 3, return; end
f = 1;
for it = 1:3
  kind = u(it).kind; if strcmp(kind, 'liter'), kind = 'litre'; end
  if ~isfield(expected, kind) || double(u(it).exponent) ~= expected.(kind), f = nan; return; end
  scale = double(u(it).multiplier) * 10^double(u(it).scale);
  if strcmp(kind, 'mole'), scale = 1000 * scale; end   % mol -> mmol
  f = f * scale ^ double(u(it).exponent);
end
