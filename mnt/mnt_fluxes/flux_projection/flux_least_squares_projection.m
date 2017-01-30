function v = flux_least_squares_projection(network,v_pre,v_sign)

% v = flux_least_squares_projection(network,v_pre)
%
% project incomplete flux vector v_pre to stationary flux space 
% v_pre flux vector with missing elements (nan) to be completed

eval(default('v_sign','[]'));  
  
[nm,nr] = size(network.N);
ind     = find(isfinite(v_pre));
 
N   = network.N(find(network.external==0),:);
dum = eye(nr); dum = dum(ind,:);
M   = full([N; dum]);
v   = inv(M'*M + 0.000001 * eye(nr)) * M' * [zeros(size(N,1),1); v_pre(ind)];
v   = project_fluxes(network.N, find(network.external), v,[],v_sign);
