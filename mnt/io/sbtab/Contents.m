% Conversion between SBtab and networks of the Metabolic Network Toolbox
%
% network_to_sbtab               - convert a network into an SBtab document
% sbtab_to_network               - convert an SBtab document into a network
% sbtab_network_state_save       - save a network with concentrations and fluxes
% sbtab_network_state_load       - load a network with concentrations and fluxes
% network_sbtab_check_for_errors - check that the Reaction and Compound tables exist
% sbtab_load_network_model_and_data - load a network model, kinetic data and state data
% sbtab_save_network_model_and_data - save a network model
%
% Metabolic states in SBtab documents (sbtab_network_state_save / _load): tables with the TableIDs
% MetaboliteConcentration, Flux, ReactionGibbsFreeEnergy and EquilibriumConstant (Quantity tables
% with columns !QuantityType, !Metabolite or !Reaction, !Value). For the corresponding SBML format,
% see 'help sbml'; for conversion, sbtab_to_sbml and sbml_to_sbtab.
