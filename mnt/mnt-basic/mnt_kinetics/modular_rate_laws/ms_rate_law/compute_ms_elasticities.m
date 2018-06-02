function E = compute_ms_elasticities(N, W, ext_ind, alpha_A, alpha_I, alpha_M, vplus, vminus, u, c)

% E = compute_ms_elasticities(N, W, ext_ind, alpha_A, alpha_I, alpha_M, vplus, vminus, u, c)

[nm,nr] = size(N);
v       = vplus  - vminus;
zeta    = vplus ./ vminus;

% -----------------------------------------------
% structure matrices and saturation values

[Mplus, Mminus, Wplus, Wminus, nm, nr, N_int] = make_structure_matrices(N,W,ext_ind);

[beta_A,gamma_A] = alpha_to_betagamma(alpha_A);
[beta_I,gamma_I] = alpha_to_betagamma(alpha_I);
[beta_M,gamma_M] = alpha_to_betagamma(alpha_M);

% -----------------------------------------------
% scaled first-order metabolite elasticities

E.sc_E_plus_c  =  Wplus .* alpha_A - Wminus .* beta_I ...
                + Mplus .* alpha_M - Mminus .* beta_M;

E.sc_E_minus_c =  Wplus .* alpha_A - Wminus .* beta_I ...
                - Mplus .* beta_M  + Mminus .* alpha_M;

E.sc_E_c       =  Wplus  .* alpha_A - Wminus .* beta_I ...
                + Mplus  .* [ alpha_M + repmat(   1./(zeta-1),1,nm)]  ...
                + Mminus .* [ alpha_M - repmat(zeta./(zeta-1),1,nm)];


% -----------------------------------------------------------------------------------
% unscaled first-order elasticities

E.un_E_plus_c   = diag(vplus)  * E.sc_E_plus_c  * diag(1./c)  ;
E.un_E_minus_c  = diag(vminus) * E.sc_E_minus_c * diag(1./c)  ;
%E.un_E_c        = E.un_E_plus_c - E.un_E_minus_c;
E.un_E_c      =  diag(v)  * [ Wplus .* alpha_A - Wminus .* beta_I ] * diag(1./c) ...
                + Mplus  .* [ diag(vplus) * alpha_M + diag(vminus) *  beta_M ] * diag(1./c) ...
                - Mminus .* [ diag(vplus) *  beta_M + diag(vminus) * alpha_M ] * diag(1./c) ;


% -----------------------------------------------
% scaled, directed second-order metabolite elasticities

E.sc_E_plus_cc   = zeros(nr,nm,nm);
E.sc_E_minus_cc  = zeros(nr,nm,nm);

diagonal_plus = Wplus .* gamma_A + Wminus .* gamma_I ...
              + Mplus .* gamma_M + Mminus .* gamma_M;

diagonal_minus = Wplus .* gamma_A + Wminus .* gamma_I ...
               + Mplus .* gamma_M + Mminus .* gamma_M;

for it = 1:nm,  
  E.sc_E_plus_cc(:,it,it)  = E.sc_E_plus_cc(:,it,it)  - diagonal_plus(:,it); 
  E.sc_E_minus_cc(:,it,it) = E.sc_E_minus_cc(:,it,it) - diagonal_minus(:,it); 
end


% -----------------------------------------------------------------------------------
% second order metabolite elasticities

E.un_E_plus_cc  = second_derivative_scaled2unscaled(E.sc_E_plus_cc,E.sc_E_plus_c,vplus,c);
E.un_E_minus_cc = second_derivative_scaled2unscaled(E.sc_E_minus_cc,E.sc_E_minus_c,vminus,c);
E.un_E_cc       = E.un_E_plus_cc - E.un_E_minus_cc;                        
E.sc_E_cc       = second_derivative_unscaled2scaled(E.un_E_cc,E.un_E_c,v,c);


% -----------------------------------------------
% first order enzyme elasticities

E.sc_E_plus_u    = speye(nr);
E.sc_E_minus_u   = speye(nr);
E.sc_E_u         = speye(nr);

E.un_E_plus_u   = diag(vplus./u);
E.un_E_minus_u  = diag(vminus./u);
E.un_E_u        = diag(v./u);


% -----------------------------------------------
% second order mixed metabolite/enzyme elasticities

E.un_E_plus_cu   = zeros(nr,nm,nr);
E.un_E_minus_cu  = zeros(nr,nm,nr);
E.un_E_cu        = zeros(nr,nm,nr);

for it = 1:nr,
  E.un_E_plus_cu(it,:,it)  = E.un_E_plus_c(it,:) / u(it);
  E.un_E_minus_cu(it,:,it) = E.un_E_minus_c(it,:)  / u(it);
  E.un_E_cu(it,:,it)       = E.un_E_c(it,:) / u(it);
end

E.sc_E_plus_cu  = zeros(nr,nm,nr);
E.sc_E_minus_cu = zeros(nr,nm,nr);
E.sc_E_cu       = zeros(nr,nm,nr);


% -----------------------------------------------
% second order enzyme elasticities

E.sc_E_plus_uu   = zeros(nr,nr,nr);
E.sc_E_minus_uu  = zeros(nr,nr,nr);
E.sc_E_plus_uu   = zeros(nr,nr,nr);
E.sc_E_minus_uu  = zeros(nr,nr,nr);
E.sc_E_uu        = zeros(nr,nr,nr);

E.un_E_plus_uu   = zeros(nr,nr,nr);
E.un_E_minus_uu  = zeros(nr,nr,nr);
E.un_E_plus_uu   = zeros(nr,nr,nr);
E.un_E_minus_uu  = zeros(nr,nr,nr);
E.un_E_uu        = zeros(nr,nr,nr);

