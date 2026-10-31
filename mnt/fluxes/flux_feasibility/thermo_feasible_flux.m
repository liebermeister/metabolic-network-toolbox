function [feasible, dmu] = thermo_feasible_flux(v, Ntot, dmu, rho, epsilon)

% [feasible, dmu] = THERMO_FEASIBLE_FLUX(v,Ntot,dmu,rho,epsilon)
%
% Test a flux vector v for EBA feasibility
%
% Criterion: existence of a vector dmu satisfying
%
%    K' * dmu = 0     (Wegscheider condition)
%    diag(sign(v)) * dmu < -rho  (Flux driven by thermodynamic force)
%
% where K is the right kernel matrix of the stoichiometric matrix Ntot.
%
% The function tests a given vector dmu for feasibility.
% If no dmu is provided as an argument, linear programming
% is used to check if such a vector exists
%
% Arguments:
%  v:       vector of fluxes (or a matrix with several flux distributions in different rows)
%  Ntot:    stoichiometric matrix, refers to all (not only to balanced) metabolites.
%  dmu:     vector of chemical potential differences (or a matrix, according to argument 'v')
%  rho:     minimal reaction Gibbs free energy to drive a non-zero flux
%  epsilon: numerical tolerance for Wegscheider condition
%
% The script requires the matlab function 'linprog' from the Optimization toolbox
%
% W. Liebermeister 2012

% Demo example:
% Ntot = [-1 1; 1 -1]; v = [1 1]';  [feasible, dmu] = thermo_feasible_flux(v,Ntot,[],0.1)
% Ntot = [-1 1; 1 -1]; v = [1 0]';  [feasible, dmu] = thermo_feasible_flux(v,Ntot,[],0.1)

if ~exist('dmu','var'),     dmu     = []; end
if ~exist('rho','var'),     rho     = 0.1; end
if ~exist('epsilon','var'), epsilon = 10^-5; end


% --------------------------------------------------------------
% If several flux vectors are given (in matrix):
% run the function separately for each of them

if size(v,2)>1
  for it = 1:size(v,2)
    if length(dmu)
      [feasible(1,it), dmu(:,it)] = thermo_feasible_flux(v(:,it),Ntot,dmu(:,it),rho,epsilon);
    else
      feasible(1,it) = thermo_feasible_flux(v(:,it),Ntot,[],rho,epsilon);
    end
  end
  return
end


% --------------------------------------------------------------
% Initialisation

nr       = size(Ntot,2);
Ktot     = null(full(Ntot)); % compute right kernel matrix
feasible = 1;

if length(dmu)
  %% dmu vector given as input: just check if v and dmu are consistent

  if max(abs(Ktot'*dmu)) > epsilon,             feasible = 0; end
  if sum( [v~=0] .* [ sign(v).*dmu > -rho] ),   feasible = 0; end

else
  %% vector dmu unknown: try to find feasible dmu vector by linear programming

  c       = ones(nr,1);
  h       = [rho * ones(nr,1); -ones(nr,1); -ones(nr,1);];    h(find(v==0)) = -1;
  G       = [-diag(sign(v));         eye(nr);      -eye(nr);];
  A       = Ktot';
  b       = zeros(size(Ktot,2),1);

  opt         = optimset('Display','off');
  opt.MaxIter = 10^10;
  [dmu, dum, exitflag] = linprog(-c, -G, -h, A, b, [], [], opt);

  if isempty(dmu) + [exitflag==-2]
    feasible = 0;
    dmu      = [];
    display(sprintf('Flux distribution is thermodynamically unfeasible at rho = %f.',rho));
  end

end
