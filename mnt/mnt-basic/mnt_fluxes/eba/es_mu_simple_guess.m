function [A, mu, exitflag] = es_mu_simple_guess(N,v,es_constraints,zv,flag_reduce_network)

% [A, mu, exitflag] = es_mu_simple_guess(N,v,es_constraints,zv,flag_reduce_network)
%
% Determine feasible reaction affinities A and chemical potentials mu from fluxes v and es_constraints

eval(default('flag_reduce_network','0','es_constraints','[]','zv','ones(size(v))'));

% Conditions:  (VERSION FOR MU):
% Thermodynamic flux direction
%   (ind_forward and ind_reverse are determined from flux directions)
%   - network.N(:,ind_forward)' * mu >= 0
%   - network.N(:,ind_reverse)' * mu <= 0
% Lower and upper bounds on mu differences
%   network.N' * mu >= dmu_min 
%   network.N' * mu <= dmu_max 
% Lower and upper bounds on mu
%   mu >= mu_min 
%   mu <= mu_max 
% Optimality criterion
%   linear weight vector could, e.g.,  be chosen as for fba
%   max = zv .* network.N' * mu

v(abs(v) < 10^-10 * max(abs(v)) ) = 0;

N_orig = N;
v_orig = v; 

if flag_reduce_network,
  %consider only subnetwork of reactions with finite fluxes
  %(some thermodynamic es_constraints can get lost!!!)
  [nm_tot, nr_tot] = size(N);
  ind_reactions   = find(v~=0);
  ind_metabolites = find(sum(abs(N(:,ind_reactions)),2));
  N = N(ind_metabolites,ind_reactions);
  v = v(ind_reactions);
  es_constraints.mu_min  = es_constraints.mu_min(ind_metabolites);
  es_constraints.mu_max  = es_constraints.mu_max(ind_metabolites);
  es_constraints.dmu_min = es_constraints.dmu_min(ind_reactions);
  es_constraints.dmu_max = es_constraints.dmu_max(ind_reactions);
  zv = zv(ind_reactions);
end

ind_forward = find(v>0);
ind_reverse = find(v<0);

M = [...
     N(:,ind_forward)'; ...
    -N(:,ind_reverse)'; ...
    -N'; ...
     N'; ...
];

isforward = [v>=0];
isreverse = [v<0];
no_dmu_max = ~isfinite(es_constraints.dmu_max);
no_dmu_min = ~isfinite(es_constraints.dmu_min);  

es_constraints.dmu_min(find([isforward .* no_dmu_min])) = -es_constraints.dmu_limit;
es_constraints.dmu_max(find([isforward .* no_dmu_max])) = -es_constraints.dmu_limit_min;
es_constraints.dmu_min(find([isreverse .* no_dmu_min])) = es_constraints.dmu_limit_min;
es_constraints.dmu_max(find([isreverse .* no_dmu_max])) = es_constraints.dmu_limit;

es_constraints.dmu_min(find([isforward .* [es_constraints.dmu_min<-es_constraints.dmu_limit]]))     = -es_constraints.dmu_limit;
es_constraints.dmu_max(find([isforward .* [es_constraints.dmu_max>-es_constraints.dmu_limit_min]])) = -es_constraints.dmu_limit_min;
es_constraints.dmu_min(find([isreverse .* [es_constraints.dmu_min<es_constraints.dmu_limit_min]]))  = es_constraints.dmu_limit_min;
es_constraints.dmu_max(find([isreverse .* [es_constraints.dmu_max>es_constraints.dmu_limit]]))      = es_constraints.dmu_limit;

b = [...
     zeros(length(ind_forward),1); ...
     zeros(length(ind_reverse),1); ...
   - es_constraints.dmu_min ; ...
     es_constraints.dmu_max ; ...
    ];

z = -N*zv;

[my_mu, val, exitflag] = lp_solve(-z,M,b,[],[],es_constraints.mu_min,es_constraints.mu_max);

if exitflag ~=1,
  exitflag
  error('No solution found');
end

%violation_ratio = sum(M*mu>b)/length(b);
my_A   = -N'*my_mu;

if flag_reduce_network,
  mu = zeros(nm_tot,1); mu(ind_metabolites) = my_mu;
  A  = -N_orig'*mu;
else,
  A  = my_A;
  mu = my_mu;
end

if [norm(N_orig' * mu + A) == 0] * isempty(find(sign(v_orig .* A)<0)),
  display('Feasible solution found');
end