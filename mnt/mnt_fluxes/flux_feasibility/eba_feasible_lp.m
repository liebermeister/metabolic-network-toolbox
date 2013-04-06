function [eba_feasible, dmu] = eba_feasible_lp(v,Ntot,dmu_fix,rho,epsilon,dmu_limit)

% [eba_feasible, dmu] = eba_feasible(v,Ntot,dmu_fix,rho,epsilon,dmu_limit)
%
% Test a flux vector v for EBA feasibility 
%
% criterion: there must exist a vector dmu such that
%         K' * dmu       = 0     (Wegscheider condition)
%    diag(sign(v)) * dmu < -rho  (Flux driven by thermodynamic force)
%               abs(dmu) <= dmu_limit 
%
% if a vector dmu_fix is given -> test this vector for feasibility
% otherwise, check by linear programming is such a vector exists
%
% rho: minimal reaction gibbs free energy to drive a non-zero flux
% epsilon: numerical tolerance for Wegscheider condition

% test: 
% Ntot = [-1 1; 1 -1]; v = [1 1]';  [eba_feasible, dmu] = eba_feasible_lp(v,Ntot,[],0.1)

eval(default('dmu_fix','[]','rho','[]','epsilon','[]','dmu_limit','1'));

if isempty(rho), rho = 10^-3; end 
if isempty(epsilon), epsilon = 10^-5; end 

% --------------------------------------------------------------
% handle several flux vectors (in matrix) 

if size(v,2)>1,
  for it = 1:size(v,2),
    if length(dmu_fix),
      [eba_feasible(1,it), dmu(:,it)] = eba_feasible_lp(v(:,it),Ntot,dmu_fix(:,it),rho,epsilon);
    else
      eba_feasible(1,it) = eba_feasible_lp(v(:,it),Ntot,[],rho,epsilon);
    end
  end
 return
end


% --------------------------------------------------------------
% initialisation

nr       = size(Ntot,2);
Ktot     = analyse_N(Ntot);
eba_feasible = 1;

if sum(isfinite(dmu_fix)) == length(v),
  %% complete dmu_fix vector given as input: just check if it works
  
  if max(abs(Ktot'*dmu_fix)) > epsilon,             eba_feasible = 0; end
  if sum( [v~=0] .* [ sign(v).*dmu_fix > -rho] ),   eba_feasible = 0; end
  dmu = dmu_fix;
  
else,
  %% try to find feasible dmu vector

  c       = ones(nr,1);
  h       = [rho * ones(nr,1);  ];    h(find(v==0)) = -1;
  G       = [-diag(sign(v));     ];
  A       = Ktot';
  b       = zeros(size(Ktot,2),1);
  if length(dmu_fix) == length(v),
    % some dmu values predefined?
    ind_finite = find(isfinite(dmu_fix));
    dd    = eye(length(dmu_fix));
    A     = [A; dd(ind_finite,:)];
    b     = [b; dmu_fix(ind_finite) ];
  end
  lb = - dmu_limit * ones(nr,1); 
  ub = dmu_limit * ones(nr,1);
  
  %%  dmu = lp236a(-c, -G, -h, A, b);

   opt         = optimset('linprog');
   opt.MaxIter = 10^10;
   opt.Display = 'off';
   [dmu, dum, exitflag] = linprog(-c, -G, -h, A, b, lb, ub, [], opt);

   if exitflag ~=1,
     exitflag
     warning('No dmu vector found');
   end
   
   if isempty(dmu) + [exitflag==-2],
     eba_feasible = 0;  display(sprintf('Flux distribution is thermodynamically unfeasible at rho = %f.',rho));
     dmu = [];
   end
  
end
