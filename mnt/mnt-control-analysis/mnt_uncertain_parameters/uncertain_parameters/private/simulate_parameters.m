% [kinetics,distr_parameters] = simulate_parameters(network, distr_parameters);
%
% draw a network's kinetic parameters from random distributions
%
% the type of kinetics is determined by distr_parameters.kinetics_type
% {'mass-action','standard'}
%
% external metabolites requires (alternatively)
%   'mu_Sext', 'sigma_Sext'   (distribution parameters for log(Sext)
%   'Sext'                     (vector of given (non-log) values, overrides distribution)
%
% formation enthalpies requires (alternatively)
%   'mu_H', 'sigma_H'           (distribution parameters for H)
%   'H'                         (vector of given (non-log) H values, overrides distribution)
%
% For mass-action kinetics:
%  rate constants (only for mass action kinetics)
%   (1) forward/backward independent:   'mu_k', 'sigma_k'
%       (distribution parameters for log rate constants)
%   (2) use equilibrium constants: sigma_prod, beta (and H computed before)
%
% For standard kinetics
% Michaelis-Menten parameters 'mu_Vm','sigma_Vm','mu_Km','sigma_Km' distribution parameters

function [kinetics,dd] = simulate_parameters(network, d);

if isfield(network,'split_reactions'),
  split_reactions = network.split_reactions;
else,
  split_reactions = 0;
end

dd = d;
if length(d.logSext_mean),
  dd.Sext = exp(d.logSext_mean + d.logSext_cov_sqrt * randn(length(find(network.external)),1));
end

switch d.kinetics_type,

  case 'mass-action',

    kinetics.type = 'mass-action';
%    kinetics.used_fwd = network.kinetics.used_fwd;
%    kinetics.used_bwd = network.kinetics.used_bwd;

    %if isfield(d,'mu_H'),
%  dd.H =  d.mu_H + d.sigma_H * randn(length(network.metabolites),1);
%end
    k = exp(d.logk_mean + d.logk_cov_sqrt * randn(size(d.logk_mean)));

    %if split_reactions,
    %  kinetics.k_fwd = k;
    %  kinetics.k_bwd = zeros(size(k));
    %  kinetics.k_fwd(find(network.used==0))=0;
    %else,
    kinetics.k_fwd = k(1:length(k)/2);
    kinetics.k_bwd = k(length(k)/2+1:end);
    %end

    if isfield(kinetics,'used_fwd'),
      kinetics.k_fwd(find(kinetics.used_fwd==0)) = 0;
      kinetics.k_bwd(find(kinetics.used_bwd==0)) = 0;
    end
%    if isfield(d,'sigma_r'),
%       [kinetics.k_fwd, kinetics.k_bwd] = simulate_k_fixed_H(network,d.H,d.beta,d.sigma_r);
%    else
%      kinetics.k_fwd = draw_from_lognormal(length(network.actions),1,d.mu_k,d.sigma_k);
%       kinetics.k_bwd = draw_from_lognormal(length(network.actions),1,d.mu_k,d.sigma_k);
%    end
%    kinetics.k_bwd(find(network.reversible==0))=0;

  case {'convenience'},
    fprintf('Warning: kinetics type not supported\n');

  case {'standard','multiple_michaelis_menten','from_network'}

   switch network.kinetics.type,
     case 'standard',
%    kinetics  = simulate_standard_kinetics(network, d,d.H);
     pars = exp(d.logk_mean + d.logk_cov_sqrt * randn(size(d.logk_mean)));
     [kinetics,Sext] = par2standard(network,[pars;dd.Sext],length(dd.Sext));
     otherwise,
       vector = exp(d.logk_mean + d.logk_cov_sqrt * randn(size(d.logk_mean)));
       kinetics = vector2parameters(network.kinetics,vector);
   end
     otherwise fprintf('Warning: kinetics type not supported\n');

end

% -----------------------------------------------------------

function [k_fwd, k_bwd, seed,H, q, products] = simulate_k_fixed_H(network,H,beta,sigma_prod,seed)

if exist('seed','var'),     randn('state',seed); end

%seed = randn('state');

% q = k_fwd / k_bwd

% r = k_fwd * k_bwd

if isfield(network,'split_reactions'), if network.split_reactions,
  q = exp(-[ network.N(:,1:length(network.actions)/2)' - network.N(:,length(network.actions)/2+1:end)']*beta*H);

    if exist('seed','var'),   r = draw_from_lognormal(length(network.actions)/2,1,0,sigma_prod,seed);
    else,                     r = draw_from_lognormal(length(network.actions)/2,1,0,sigma_prod);
    end
    k_fwd = [sqrt(r.*q);sqrt(r./q)];
    k_bwd = 0*k_fwd;
  end
else,

q = exp(-network.N'*beta*H);

if exist('seed','var'),   r = draw_from_lognormal(length(network.actions),1,0,sigma_prod,seed);
else,                     r = draw_from_lognormal(length(network.actions),1,0,sigma_prod);
end

  k_fwd = sqrt(r.*q);
  k_bwd = sqrt(r./q);
end

% ------------------------------------------------------------

function [m,seed] = draw_from_lognormal(ni,nk,mean,sigma,seed)

if exist('seed','var'), randn('state',seed); end

m    = exp(mean+sigma*randn(ni,nk));
seed = randn('state');
