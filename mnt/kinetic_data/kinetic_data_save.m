function M = kinetic_data_save(kinetic_data, network, filename, options_style, verbose)

% M = kinetic_data_save(kinetic_data,network,filename,options_style,verbose)
%
% Save a kinetic data struct to an SBtab file (or display it, if no filename is given)
%
% also see kinetic_data_load, kinetic_data_to_sbtab (conversion into an SBtab table object)

eval(default('filename','[]','options_style','[]','verbose','0'));

if isempty(options_style)
  options_style=0;
end

[~, M] = kinetic_data_to_sbtab(kinetic_data, network);

if length(filename)
  if verbose
    display(sprintf('Writing kinetic data to SBtab file %s ',filename));
  end
  mytable(M,options_style,filename);
else
  mytable(M,options_style)
end
