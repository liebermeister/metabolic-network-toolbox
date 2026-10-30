function sbtab_save_network_model_and_data(my_object, filename)

% sbtab_save_network_model_and_data(my_object, filename)
%
% Save a network model (field my_object.network) to an SBtab file
%
% Saving state data and options is not implemented yet.
% For loading network model (including kinetic data) and state data, see: sbtab_load_network_model_and_data

if ~isfield(my_object,'network')
  error('The object has no field "network"');
end

s = network_to_sbtab(my_object.network);

sbtab_document_save_to_one(s,filename);
