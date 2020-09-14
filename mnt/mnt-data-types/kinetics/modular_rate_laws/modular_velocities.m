function [v, v_plus, v_minus, D, regulation_term] = modular_velocities(kinetic_law, N, W, ind_ext, u, c, KA, KI, KM, KV, Keq, h, Mplus, Mminus, Wplus, Wminus, nm, nr)

% [v, v_plus, v_minus] = modular_velocities(kinetic_law,N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h,Mplus, Mminus, Wplus, Wminus, nm, nr)

% kinetic_law  string:  {'cs','ms','ds','rp','fd'}, corresponding to rate laws in modular rate law paper
% nm       # metabolites
% nr       # reactions
% N        mn x nr Stoichiometric matric
% W        nr x nm Allosteric regulation matrix (+1 for activator, -1 for inhibitor) 
% ind_ext  vector of indices of external metabolites (numbers between 1 and nm)
% u        nr x 1 column vector of enzyme concentrations
% c        nm x 1 column vector of enzyme concentrations
% KA       nr x nm sparse matrix of KA values (activation)
% KI       nr x nm sparse matrix of KI values (inhibition)
% KM       nr x nm sparse matrix of KM values (Michaelis-Menten)
% KV       nr x 1 column vector of KV values (velocity constants, defined as sqrt(Kcat_forward * Kcat_reverse)
% Keq      nr x 1 column vector of equilibrium constants
% h        nr x 1 column vector of "cooperativity factors" (usually, they can all be set to 1)
% Mplus    nr x nm matrix of forward molecularities (ie |N|, but only for substrate elements)
% Mminus   nr x nm matrix of reverse molecularities (ie N, but only for product elements)
% Wplus    nr x nm matrix of activations (ie W, but only for activation elements)
% Wminus   nr x nm matrix of inhibitions (ie |W|, but only for inhibition elements)
% 
% Note that Mplus, Mminus, Wplus, Wminus, nm, nr can be computed from N,W,ind_ext,h by using the matlab function make_structure_matrices
    
c(find(c<10^-14)) = 10^-14;

%% Global variables to speed up this function; see ecm_one_run.m and convex_parameter_estimation.m
global global_structure_matrices 

if exist('global_structure_matrices','var'),
  global Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm
end

if isempty(Mplus),
  [Mplus, Mminus, Wplus, Wminus, nm, nr, N_int,ind_M,ind_Wp,ind_Wm] = make_structure_matrices(N,W,ind_ext,h);
end

[log_alpha_A, log_beta_A] = k_to_log_alpha(KA,c);
[log_alpha_I            ] = k_to_log_alpha(KI,c);
[Kplus,Kminus]            = ms_compute_Kcat(N,KM,KV,Keq);

log_c_by_k        = repmat(log(c)',nr,1);
log_c_by_k(ind_M) = log_c_by_k(ind_M) - log(KM(ind_M)+10^-15);

theta_plus  = exp( sum( Mplus  .* log_c_by_k, 2) );
theta_minus = exp( sum( Mminus .* log_c_by_k, 2) );

switch kinetic_law,  
  case {'cs','ms'},
    psi_plus    = exp( sum( Mplus  .* log(1+exp(log_c_by_k)) , 2) );
    psi_minus   = exp( sum( Mminus .* log(1+exp(log_c_by_k)) , 2) );    
end

switch kinetic_law,  
  case 'cs', D = psi_plus + psi_minus - 1;
  case 'ms', D = psi_plus .* psi_minus;
  case 'ds', D = theta_plus + theta_minus + 1;
  case 'rp', D = 1;
  case 'fd', D = sqrt(theta_plus .* theta_minus);
end

% complete allosteric regulation
regulation_term = exp( sum( Wplus .* log_beta_A + Wminus .* log_alpha_I, 2) );

v_plus  = u .* regulation_term .* Kplus  .* theta_plus  ./ D;
v_minus = u .* regulation_term .* Kminus .* theta_minus ./ D;

v_plus  = real(full(v_plus));
v_minus = real(full(v_minus));
v       = v_plus - v_minus;
