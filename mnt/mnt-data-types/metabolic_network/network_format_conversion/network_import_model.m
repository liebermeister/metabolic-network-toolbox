function network = network_import_model(model_file,options);

eval(default('options','struct'));

options_default = struct('kinetic_law','cs');
options = join_struct(options_default,options);

% REDUNDANT WITH FUNCTION network_import!!  
  
switch model_file(end-3:end),
  case '.xml',
    format = 'sbml';
  case '.tsv',
    format = 'sbtab';
  otherwise
    error(sprintf('File %s: unknown file format',model_file));
end

switch format,
  case 'sbml'
    network = network_sbml_import(model_file);
  case 'sbtab'
    network = sbtab_to_network(model_file, options);
end
