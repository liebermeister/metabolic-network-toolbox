function network = network_import_model(model_file,options);

% network = network_import_model(model_file,options)
%
% Load a network model from a file (see also network_import)
%
% model_file: one of the following
%   o SBML file name (extension '.xml')
%   o SBtab file name (extension '.tsv')
%   o SBtab document object (Matlab struct with field 'tables'; tables 'Reaction' and 'Compound')
%   o network struct (Metabolic Network Toolbox; returned unchanged)

eval(default('options','struct'));

options_default = struct('kinetic_law','cs');
options = join_struct(options_default,options);

% REDUNDANT WITH FUNCTION network_import!!

if isstruct(model_file)
  if isfield(model_file,'tables')
    % SBtab document object
    network = sbtab_to_network(model_file, options);
  elseif isfield(model_file,'N') && isfield(model_file,'metabolites')
    % network struct
    network = model_file;
  else
    error('Unknown model data structure (expected an SBtab document object or a network struct)');
  end
  return
end

switch model_file(end-3:end)
  case '.xml'
    format = 'sbml';
  case '.tsv'
    format = 'sbtab';
  otherwise
    error(sprintf('File %s: unknown file format',model_file));
end

switch format
  case 'sbml'
    network = network_sbml_import(model_file);
  case 'sbtab'
    network = sbtab_to_network(model_file, options);
end
