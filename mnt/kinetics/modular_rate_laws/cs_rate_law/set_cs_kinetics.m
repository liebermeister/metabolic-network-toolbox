function kinetics = set_cs_kinetics(network,parameters)

% SET_CS_KINETICS Construct common saturable (cs) kinetics
%
% kinetics = set_cs_kinetics(network,parameters)
%
% parameters: see set_ms_kinetics
%
% See also: set_kinetics, set_ms_kinetics

eval(default('parameters','struct'));

kinetics = set_ms_kinetics(network,parameters);

kinetics.type = 'cs';
