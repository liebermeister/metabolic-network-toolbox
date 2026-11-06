% TEMPORAL_CONTROL_COEFFICIENTS Time-dependent control coefficients after a step change of enzyme activities
%
% [CS_t, CJ_t, CS_t_scaled, CJ_t_scaled] = temporal_control_coefficients(N, Ec, external, t, c, v, mu)
%
% Linear response of concentrations and fluxes, starting from a stable steady state, to a step
% change of reaction rates (enzyme activities) at time 0. In the linear approximation, with the
% Jacobian M = NR * Ec_int * L of the independent internal metabolites,
%
%   CS(t) = L * (expm(M t) - I) * inv(M) * NR       (concentrations, unscaled)
%   CJ(t) = I + Ec_int * CS(t)                       (fluxes, unscaled)
%
% For t = 0, CS = 0 and CJ = I; for large t, they approach the steady-state control coefficients
% (see control_coefficients). The coefficients describe responses to perturbations of the
% reaction rates (dv = du for unchanged concentrations), as in control_coefficients.
%
% Inputs
%  N          stoichiometric matrix (all metabolites)
%  Ec         matrix of unscaled reaction elasticities (reactions x metabolites)
%  external   binary vector marking external metabolites
%  t          vector of time points
%  c, v       (optional) concentrations and fluxes in the steady state, for the scaled coefficients
%  mu         (optional) growth rate, i.e. dilution rate of internal metabolites (default 0), as in
%             control_coefficients: with mu > 0, all internal metabolites are diluted (Jacobian
%             N_int * Ec_int - mu * I), and there are no conservation relations
%
% Outputs
%  CS_t        unscaled concentration control coefficients, array (metabolites x reactions x time points);
%              rows of external metabolites are zero
%  CJ_t        unscaled flux control coefficients, array (reactions x reactions x time points)
%  CS_t_scaled scaled coefficients d ln c / d ln u (requires c and v; external metabolites: 0)
%  CJ_t_scaled scaled coefficients d ln v / d ln u (requires v; rows of reactions with zero flux: nan)
%
% Scaled concentration coefficients can be used for small relative changes of enzyme activities:
%  d ln c(t) = CS_t_scaled(:,:,i) * d ln u
%
% See also: control_coefficients, time_response_coefficients (nonlinear trajectories, requires kinetics)

function [CS_t, CJ_t, CS_t_scaled, CJ_t_scaled] = temporal_control_coefficients(N, Ec, external, t, c, v, mu)

if ~exist('mu','var') || isempty(mu), mu = 0; end

if any(~isfinite(Ec(:))), error('Elasticity matrix contains non-finite values'); end

[nm, nr] = size(N);
ind_int  = find(external(:) == 0);
N_int    = N(ind_int,:);
Ec_int   = Ec(:,ind_int);

if mu == 0
  [L, NR] = reduce_N(N_int);
  M       = NR * Ec_int * L;
else
  % with dilution, all internal metabolites are independent
  L  = eye(length(ind_int));
  NR = N_int;
  M  = N_int * Ec_int - mu * eye(length(ind_int));
end

if rank(M) < size(M,1)
  error('The Jacobian matrix is singular');
end
if max(real(eig(M))) > 1e-10
  warning('The steady state is unstable: the coefficients grow with time');
end

M_inv_NR = M \ NR;
n_t      = length(t);
CS_t     = zeros(nm, nr, n_t);
CJ_t     = zeros(nr, nr, n_t);

for it = 1:n_t
  CS_int         = L * (expm(M * t(it)) - eye(size(M))) * M_inv_NR;
  CS_t(ind_int,:,it) = CS_int;
  CJ_t(:,:,it)   = eye(nr) + Ec_int * CS_int;
end

if nargout > 2
  if ~exist('c','var') || ~exist('v','var') || isempty(c) || isempty(v)
    error('Scaled coefficients require the concentrations c and the fluxes v');
  end
  c = c(:); v = v(:);
  CS_t_scaled = zeros(size(CS_t));
  CJ_t_scaled = zeros(size(CJ_t));
  v_inv = 1 ./ v;
  v_inv(v == 0) = nan;
  for it = 1:n_t
    CS_t_scaled(ind_int,:,it) = diag(1 ./ c(ind_int)) * CS_t(ind_int,:,it) * diag(v);
    CJ_t_scaled(:,:,it)       = diag(v_inv) * CJ_t(:,:,it) * diag(v);
  end
end
