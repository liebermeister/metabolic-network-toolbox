function sbtab_object = sbml_to_sbtab(sbml_object, filename)

% sbtab_object = sbml_to_sbtab(sbml_object, filename)
%
% Convert an SBML model (libSBML Matlab struct, or SBML file name) into an SBtab document.
% If 'filename' is given, the SBtab document is also written to a file.
%
% Metabolic state: if the SBML model contains the local reaction parameters 'reference_flux',
% 'equilibrium_constant' or 'reaction_gibbs_free_energy' (SBML format of the Metabolic Network
% Toolbox for models with a metabolic state, written by sbml_network_state_save), the model is read
% with sbml_network_state_load, and the SBtab document contains the tables MetaboliteConcentration
% (from the initial concentrations of the species), Flux, ReactionGibbsFreeEnergy and
% EquilibriumConstant (as written by sbtab_network_state_save). Otherwise, only the network is
% converted (network_sbml_import, network_to_sbtab).
%
% This function requires the Metabolic Network Toolbox and libSBML (TranslateSBML)
%
% usage example:
% sbml_object = sbtab_load_example('sbml/escherichia-coli_noor_2016.xml','sbml');
% sbtab_object = sbml_to_sbtab(sbml_object);
%
% See also: sbtab_to_sbml, sbml_network_state_load, sbtab_network_state_save

eval(default('filename','[]'));

if ischar(sbml_object) || isstring(sbml_object)
  sbml_object = TranslateSBML(char(sbml_object), 0, 0);
end

if sbml_contains_state(sbml_object)
  warning_state = warning('off', 'all');   % missing state values are fine here
  [network, c, v, dmu, keq] = sbml_network_state_load(sbml_object);
  warning(warning_state);
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
  kl = s.reaction(it).kineticLaw;
  if isempty(kl), continue; end
  if isfield(kl,'localParameter'), lp = kl.localParameter; else, lp = kl.parameter; end
  if ~isempty(lp) && any(ismember({lp.id}, {'reference_flux', 'equilibrium_constant', 'reaction_gibbs_free_energy'}))
    flag = 1;
    return
  end
end
