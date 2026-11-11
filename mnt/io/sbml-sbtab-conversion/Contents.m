% SBML <-> SBtab conversion via MNT network objects
%
% sbml_to_sbtab         - convert an SBML model into an SBtab document
% sbtab_to_sbml         - convert an SBtab model into an SBML model
%
% Both functions also convert metabolic states (concentrations, fluxes, equilibrium constants,
% reaction Gibbs free energies). In SBtab, the state is stored in the tables MetaboliteConcentration,
% Flux, EquilibriumConstant and ReactionGibbsFreeEnergy; in SBML, in the MNT format for models with a
% metabolic state: concentrations as initial concentrations of the species, and fluxes, equilibrium
% constants and Gibbs free energies as local parameters 'reference_flux', 'equilibrium_constant'
% and 'reaction_gibbs_free_energy' of each reaction (see 'help sbml' and sbml_network_state_save).
