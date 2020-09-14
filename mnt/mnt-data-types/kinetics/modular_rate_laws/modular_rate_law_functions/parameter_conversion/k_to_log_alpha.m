function [log_alpha_X,log_beta_X] = k_to_log_alpha(KX,c,hill);

%function [log_alpha_X,log_beta_X] = k_to_log_alpha(KX,c,hill);
%
%Let nm x nr be the size of the stoichiometric matrix 
%
% Input:
% KX       nr x nm sparse matrix of KM values, (Michaelis-Menten), KA values (activation), or KI values (inhibition)
% c        mn x 1 column vector of metabolite concentrations
% hill     nr x 1 column vector of "cooperativity factors" (usually, they can all be set to 1)
%
% Output:
% log_alpha_X: nr x nm sparse matrix of log saturation values, where alpha = 1 / (1 + c/KX)
% log_beta_X:  nr x nm sparse matrix of log saturation values, where beta (c/KX) / (1 + c/KX)

if ~exist('hill','var'), hill = 1; end

[nr,nm] = size(KX);
ind_KX  = find(KX);

C                     = repmat(c',nr,1);
log_alpha_X           = sparse(nr,nm); 
log_alpha_X(ind_KX)   = - log( 1 + [C((ind_KX)) ./ KX(ind_KX)].^hill  );

if nargout >1,
  log_beta_X          = sparse(nr,nm); 
  log_beta_X(ind_KX)  = log( [C((ind_KX)) ./ KX(ind_KX)].^hill  ) + log_alpha_X(ind_KX);
end

