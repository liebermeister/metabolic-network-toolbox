function alpha_A = beta_to_alpha(beta_A)

% [alpha_A] = beta_to_alpha(beta_A)

[nr,nm] = size(beta_A);

alpha_A = sparse(nr,nm);
alpha_A(find(beta_A)) = 1-beta_A(find(beta_A));
