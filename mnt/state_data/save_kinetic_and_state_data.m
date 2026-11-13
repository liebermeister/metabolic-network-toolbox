function save_kinetic_and_state_data(network, kinetic_data, state_data, filename_kinetic_data, filename_state_data, consistent, standard_concentration)

% save_kinetic_and_state_data(network, kinetic_data, state_data, filename_kinetic_data, filename_state_data, consistent, standard_concentration)
%
% Save kinetic data and state data
%
% standard_concentration (optional, default '1 mM'): standard concentration ('1 mM' or '1 M') for
% equilibrium constants and standard Gibbs free energies in the kinetic data file
%
% For saving network model (with kinetic data) and state data, see: sbtab_save_network_model_and_data
% For loading kinetic data and state data, see: sbtab_load_network_model_and_data

eval(default('consistent','0','standard_concentration','''1 mM'''));

display(sprintf('Writing a copy of kinetic data used:'));
kinetic_data_save(kinetic_data,network,filename_kinetic_data,[],0,standard_concentration);

display(sprintf('Writing a copy of state data used:'));
state_data_save(state_data,network,filename_state_data,struct('consistent',consistent));
