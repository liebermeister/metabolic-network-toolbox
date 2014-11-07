function v = pmf(network,fba_constraints,benefit,v_start)

% v = pmf(network,fba_constraints,benefit,v_start)
%
% Principle of minimal sum of absolute fluxes, solves
%
% min = sum(abs(v))
% where benefit = z' * v 
%       N_internal * v = 0  
%       v_min <= v <= v_max
%
% fba_constraints: see fba_default_options
%
%  fba_constraints.zv:       linear weights in the objective function  
%  fba_constraints.v_min:    vector of lower bounds
%  fba_constraints.v_max:    vector of upper bounds
%  fba_constraints.v_sign:   vector of signs, overrides v_min and v_max
%  fba_constraints.v_fix:    vector of fixed fluxes, overrides everything else
%  fba_constraints.ext_sign: sign vector for external metabolite production
%
% benefit: predefined value of the objective
% v_start: (optional) initial flux vector to be optimised

[nm,nr] = size(network.N);

fba_constraints = fba_update_constraints(fba_constraints);

ind_fix = find(isfinite(fba_constraints.v_fix));
v_min = fba_constraints.v_min;
v_max = fba_constraints.v_max;
K    = null(full(network.N(find(network.external==0),:)));
zv   = fba_constraints.zv;

if exist('v_start','var'),
  v0 = v_start;
else,
  v0 = 0.5 * [v_min+v_max];
  v0(ind_fix) = fba_constraints.v_fix(ind_fix);
end
v0red = pinv(full(K)) * v0;

% min = sum(abs(K*vred))
% where K * vred > v_min
%       K * vred < v_max
%       K(fixed,:) * vred = v_fix
%       zv' * K * vred = benefit;

A    = [-K; K];
b    = [-v_min; v_max];
A_eq = [K(ind_fix,:); zv' * K];
b_eq = [fba_constraints.v_fix(ind_fix); benefit];

g1   = inline('sum(abs(K*vred))','vred','K');
opt  = optimset('Display','off','MaxIter',10^10,'MaxFunEvals',10^10,'Algorithm','active-set');

% check
%  v0 >= v_min 
%  v0 <= v_max
%  v0(ind_fix) == fba_constraints.v_fix(ind_fix)
%  zv' * v0 == benefit
 
% v0red = 0 * v0red;
% A    * v0red <= b
% A_eq * v0red == b_eq

[vred,fval,exitflag] = fmincon(@(vred) g1(vred,K), v0red, A, b, A_eq, b_eq, [],[],[],opt);

if exitflag<=0,
  exitflag
  warning('No solution found in FBA with minimal fluxes'); 
end

v = K * vred;
