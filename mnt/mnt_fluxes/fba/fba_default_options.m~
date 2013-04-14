function fba_constraints = fba_default_options(network)

% fba_constraints = fba_default_options(network);

[nm,nr] = size(network.N);

fba_constraints.zv             = nan*ones(nr,1);
fba_constraints.v_min          = -ones(nr,1);
fba_constraints.v_max          = ones(nr,1);
fba_constraints.v_fix          = nan*ones(nr,1);
fba_constraints.v_sign         = nan*ones(nr,1);
fba_constraints.ext_sign       = nan*ones(sum(network.external),1);
fba_constraints.concentrations = ones(nm,1);