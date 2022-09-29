function network = network_import(network_file,options)
  
% network = network_import(network_file)
% convenience function for loading a network from either an SBML or SBtab file (detected by file extension)  
  
eval(default('options','struct'));
  
switch network_file(end-3:end),
  case {'.tsv','.csv'},
    network = sbtab_to_network(network_file,options);
  case '.xml',
    network = network_sbml_import(network_file);
  otherwise 
    error('Unexpected file extension');
end
