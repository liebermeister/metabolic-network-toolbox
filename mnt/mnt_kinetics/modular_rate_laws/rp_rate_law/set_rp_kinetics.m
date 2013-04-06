function kinetics = set_rp_kinetics(network,parameters)

% kinetics = set_rp_kinetics(network,parameters)

if exist('parameters','var'),
  kinetics = set_ms_kinetics(network,parameters);
else,
  kinetics = set_ms_kinetics(network);
end

kinetics.type = 'rp';