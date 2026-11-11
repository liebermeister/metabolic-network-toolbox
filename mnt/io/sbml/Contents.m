% Import and export of SBML models (requires the Matlab interface of libSBML)
%
% network_sbml_import                  - import a network from SBML
% network_sbml_export                  - export a network to SBML
% network2sbml                         - convert a network to an SBML structure
% network_adjust_names_for_sbml_export - adjust names to the SBML syntax
% get_formula                          - rate law formula of a reaction
% sbml_network_state_save              - save a network and a metabolic state (c, v, dmu, keq) to SBML
% sbml_network_state_load              - load a network and a metabolic state (c, v, dmu, keq) from SBML
%
% SBML format for models with a metabolic state (sbml_network_state_save / _load)
%
% SBML has no elements for flux data or equilibrium constants. MNT therefore uses standard SBML
% Level 3 Version 2 elements with fixed ids, which any SBML tool can read:
%   - concentrations: initialConcentration of each <species> (default unit mM)
%   - fluxes, equilibrium constants, reaction Gibbs free energies: <localParameter> elements in
%     the <kineticLaw> of each <reaction>, with the ids
%        reference_flux              SBO:0000612 (rate of reaction), default unit mM/s
%        equilibrium_constant        SBO:0000281, unit mM^(sum of stoichiometric coefficients)
%        reaction_gibbs_free_energy  SBO:0000617 (Gibbs free energy of reaction), unit kJ/mol
%     Each parameter has a unit definition and a note explaining its meaning; the model notes
%     describe the format. All values refer to the direction of the reaction formula.
%   - kinetic laws: for networks without kinetics, the constant rate
%     reference_flux * flux_reference_volume (global parameter, 1 litre), so that a simulation
%     stays in the (stationary) reference state; existing kinetics are written as usual, and the
%     three parameters are added to the kinetic law without being used in the rate law
%   - compartments (network.compartments, compartment_sizes, metabolite_compartments), database
%     identifiers (KEGG, BiGG, ChEBI) as identifiers.org annotations, and the original metabolite
%     and reaction ids (name attributes; ids that are no valid SBML ids are adjusted)
% The concentration and time units can be chosen (options concentration_unit, time_unit).
