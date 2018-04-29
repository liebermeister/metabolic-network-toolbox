function [v, v_plus, v_minus, D, regulation_term ] = modular_velocities(kinetic_law,N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h,Mplus, Mminus, Wplus, Wminus, nm, nr)

% [v, v_plus, v_minus] = modular_velocities(kinetic_law,N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h,Mplus, Mminus, Wplus, Wminus, nm, nr)

c(find(c<10^-14)) = 10^-14;

if ~exist('Mplus','var'),
  [Mplus, Mminus, Wplus, Wminus, nm, nr] = make_structure_matrices(N,W,ind_ext,h);
end

ind_N  = find(Mplus + Mminus);
ind_Wp = find(Wplus);
ind_Wm = find(Wminus);
 
[log_alpha_A, log_beta_A] = k_to_log_alpha(KA,c);
[log_alpha_I            ] = k_to_log_alpha(KI,c);
[Kplus,Kminus]            = ms_compute_Kcat(N,KM,KV,Keq);

log_c_by_k        = repmat(log(c)',nr,1);
log_c_by_k(ind_N) = log_c_by_k(ind_N) - log(KM(ind_N)+10^-15);

theta_plus  = exp( sum( Mplus  .* log_c_by_k, 2) );
theta_minus = exp( sum( Mminus .* log_c_by_k, 2) );

switch kinetic_law,  
  case {'cs','ms'},
    psi_plus    = exp( sum( Mplus  .* log(1+exp(log_c_by_k)) , 2) );
    psi_minus   = exp( sum( Mminus .* log(1+exp(log_c_by_k)) , 2) );    
end

% complete allosteric regulation
regulation_term = exp( sum( Wplus .* log_beta_A + Wminus .* log_alpha_I, 2) );

switch kinetic_law,  
  case 'cs', D = psi_plus + psi_minus - 1;
  case 'ms', D = psi_plus .* psi_minus;
  case 'ds', D = theta_plus + theta_minus + 1;
  case 'rp', D = 1;
  case 'fd', D = sqrt(theta_plus .* theta_minus);
end

v_plus  = u .* regulation_term .* Kplus  .* theta_plus  ./ D;
v_minus = u .* regulation_term .* Kminus .* theta_minus ./ D;

v_plus  = real(full(v_plus));
v_minus = real(full(v_minus));
v       = v_plus - v_minus;

%D
%[Kplus .* theta_plus, Kminus .* theta_minus]
%[v_plus, v_minus]