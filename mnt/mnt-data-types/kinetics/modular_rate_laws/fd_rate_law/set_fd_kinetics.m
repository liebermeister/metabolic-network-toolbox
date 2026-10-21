function kinetics = set_fd_kinetics(network,parameters)

if exist('parameters','var'),
  kinetics = set_ms_kinetics(network,parameters);
else,
  kinetics = set_ms_kinetics(network);
end

kinetics.type = 'fd';