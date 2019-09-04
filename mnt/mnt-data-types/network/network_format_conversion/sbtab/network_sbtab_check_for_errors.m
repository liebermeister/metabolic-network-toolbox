function network_sbtab_check_for_errors(my_sbtab)

error = '';

if ~isfield(my_sbtab.tables,'Reaction'), error('Reaction table missing'); end
if ~isfield(my_sbtab.tables,'Compound'), error('Compound table missing'); end
