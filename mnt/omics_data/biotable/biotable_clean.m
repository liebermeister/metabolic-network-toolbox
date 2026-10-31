function data_table = biotable_clean(data_table)

% data_table = biotable_clean(data_table)
%
% Remove the fields 'metainfo' and 'file' from a biotable

if isfield(data_table,'metainfo'), data_table = rmfield(data_table,'metainfo'); end
if isfield(data_table,'file'),     data_table = rmfield(data_table,'file'); end
