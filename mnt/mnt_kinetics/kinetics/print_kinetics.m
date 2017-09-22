% print_kinetics(kinetics,verbose)

function print_kinetics(kinetics,verbose)

if ~exist('verbose'), verbose =0; end

fprintf('Type: %s\n\n',kinetics.type);

switch kinetics.type,
  
  case 'numeric',
  
    for it=1:length(kinetics.reactions),
      fprintf('Reaction %d:\t%s  \n',it,strrep(kinetics.reactions{it}.string,' ',''));
if verbose,
      fprintf('(with ');
      for it2=1:length(kinetics.reactions{it}.parameters),
	fprintf(' %s = %f ',kinetics.reactions{it}.parameters{it2}.name,kinetics.parameter_values(kinetics.reactions{it}.parameters{it2}.index));
      end
      fprintf(')\n');
end
    end
    fprintf('\n');
    for it=1:length(kinetics.parameters),
      fprintf('Parameter %s\t= %f\n',kinetics.parameters{it},kinetics.parameter_values(it));
    end

  case 'standard',
      for it=1:length(kinetics.reactions),
	  fprintf('Reaction %d: %s',it,kinetics.reactions{it}.type);
	  for itt=1:length(kinetics.reactions{it}.parameters),
	    fprintf(' %s=%f',kinetics.reactions{it}.parameters{itt},getfield(kinetics.reactions{it},kinetics.reactions{it}.parameters{itt}));
	  end
	  fprintf('\n');
      end
      
  otherwise fprintf('Warning: kinetics type not supported\n');

end

