function [dmu, success_flag] = sample_feasible_dmu(N,ind_ext,v,stm_constraints,stm_options)

% [mu, success_flag] = sample_feasible_mu(N, ind_ext, v, stm_constraints, stm_options)
%
% Sample dmu vectors that agree with a given flux vector v
% and with further stm_constraints given in 'stm_constraints':
%
% Input: dmu_min, dmu_max, sign(v), rho, dmu_limit, dmu_fix, and dmu_eqconstraint
%
% But DISREGARD all other entries in 'stm_constraints' (about mu, c, c0, Keq)
%
% Compute some dmu vectors at the boundaries of the feasible region
% set by the (previously sampled or determined) flux directions
% compute the corresponding mu vectors using the pseudoinverse of N'

dmu          = [];
success_flag = 1;
[nm,nr]      = size(N);

if ~isfield(stm_options,'seed'), stm_options.seed = 0; end
if ~isnan(stm_options.seed), rng(stm_options.seed); end

epsilon     = RT/stm_constraints.rho;
v_signs     = sign(v);

% update lower and upper bounds for dmu

dmu_min = stm_constraints.dmu_min;
dmu_max = stm_constraints.dmu_max;
dmu_min(~isfinite(dmu_min))                   = - stm_constraints.dmu_limit;
dmu_max(~isfinite(dmu_max))                   =   stm_constraints.dmu_limit;
dmu_max(v_signs>0)                            = - epsilon;
dmu_min(v_signs<0)                            =   epsilon;
dmu_min(dmu_min < - stm_constraints.dmu_limit) = - stm_constraints.dmu_limit;
dmu_max(dmu_max >   stm_constraints.dmu_limit) =   stm_constraints.dmu_limit;

ee   = eye(length(v));

A_ineq = [  ee; ...
          - ee; ...
         ];

b_ineq = [ dmu_max; ...
          -dmu_min; ...
         ];


if isempty(stm_constraints.dmu_eqconstraint)
  stm_constraints.dmu_eqconstraint = struct('matrix',zeros(0,nr),'vector',zeros(0,1));
end

A_eq = [ee(find(isfinite(stm_constraints.dmu_fix)),:); ...
        stm_constraints.dmu_eqconstraint.matrix;...
       ];

b_eq = [ stm_constraints.dmu_fix(find(isfinite(stm_constraints.dmu_fix))); ...
         stm_constraints.dmu_eqconstraint.vector; ...
       ];

[dmu, ~, ~, success_flag] = convex_sampling(A_ineq, b_ineq, A_eq, b_eq, stm_options.n_dmu_samples);

% check inequality constraints and equality constraints
% if success_flag,
%   [A_ineq * dmu, b_ineq]
%   [A_eq * dmu, b_eq]
% end
