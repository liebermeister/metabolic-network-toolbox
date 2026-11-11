% Import and export of metabolic network models
%
%   sbml                    - SBML files; requires the Matlab interface of libSBML
%   sbtab                   - SBtab files; requires SBtab for Matlab
%   sbml-sbtab-conversion   - conversion between SBML and SBtab
%   cobra                   - COBRA models
%
% network_import        - load a network from a file
% network_import_model  - load a network with kinetics from a file
%
% Models with a metabolic state (concentrations c, fluxes v, equilibrium constants keq, reaction
% Gibbs free energies dmu) can be saved and loaded in both formats:
%   SBtab: sbtab_network_state_save, sbtab_network_state_load
%          (tables MetaboliteConcentration, Flux, EquilibriumConstant, ReactionGibbsFreeEnergy)
%   SBML:  sbml_network_state_save, sbml_network_state_load
%          (MNT's SBML format for models with a metabolic state, SBML Level 3 Version 2:
%           concentrations = initial concentrations of the species; fluxes, equilibrium constants
%           and Gibbs free energies = local parameters 'reference_flux', 'equilibrium_constant',
%           'reaction_gibbs_free_energy' in the kinetic law of each reaction, with units, SBO terms
%           and notes; see 'help sbml')
% sbtab_to_sbml and sbml_to_sbtab convert models with states between the two formats.
% mnt_load_model_and_state (folder model_and_state) loads a model with its state from a .mat,
% SBtab or SBML file.
