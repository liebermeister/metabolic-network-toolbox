% Models with a reference state (network struct, fluxes v, concentrations c, equilibrium constants keq)
%
% Loading
%   mnt_load_model_and_state    - Load model and reference state from a .mat, SBtab or SBML file (type from the extension)
%
% Data mapping
%   mnt_map_state_data          - Read flux, concentration and enzyme data from SBtab files and map them to a model
%   mnt_keq_from_kinetic_data   - Equilibrium constants for a model by parameter balancing of kinetic/thermodynamic data
%
% Both data mapping functions convert the data to the reaction orientation of the model first
% (sbtab_orient_data_files, folder data_orientation) unless options.orient_data = 0.
% mnt_keq_from_kinetic_data calls parameter_balancing (parameter balancing toolbox).
% Reference states are constructed from the mapped data by stm_construct_reference_state (STM toolbox).
