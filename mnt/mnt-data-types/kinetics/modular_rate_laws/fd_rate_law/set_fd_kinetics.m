function kinetics = set_fd_kinetics(network,parameters)

% kinetics = set_fd_kinetics(network,parameters)
%
% Construct a kinetics structure for the force-dependent (fd) modular rate law; see set_ms_kinetics

if exist('parameters','var')
  kinetics = set_ms_kinetics(network,parameters);
else
  kinetics = set_ms_kinetics(network);
end

kinetics.type = 'fd';
