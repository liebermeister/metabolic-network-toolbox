function save_kinetic_and_state_data(network, kinetic_data, state_data, filename_kinetic_data, filename_state_data,consistent);

% save_kinetic_and_state_data(kinetic_data, state_data, filename_kinetic_data, filename_state_data);
%
% Save kinetic data and state data
%
% For saving network model (with kinetic data) and state data, see: sbtab_save_network_model_and_data
% For loading kinetic data and state data, see: sbtab_load_network_model_and_data

eval(default('consistent','0'));

display(sprintf('Writing a copy of kinetic data used:'));
kinetic_data_save(kinetic_data,network,filename_kinetic_data);

display(sprintf('Writing a copy of state data used:'));
state_data_save(state_data,network,filename_state_data,struct('consistent',consistent));
