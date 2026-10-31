function errors = network_sbtab_check_for_errors(my_sbtab)

% errors = network_sbtab_check_for_errors(my_sbtab)
%
% Check that an SBtab document contains the tables needed for a network model
% (an error is raised if the Reaction or Compound table is missing)

errors = '';

if ~isfield(my_sbtab.tables,'Reaction'), error('Reaction table missing'); end
if ~isfield(my_sbtab.tables,'Compound'), error('Compound table missing'); end
