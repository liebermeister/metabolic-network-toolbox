function [this_mu,fval,e] = fminconIPasQPsolver(H,f,A,b,mu_lower,mu_upper,mu_init)

% [this_mu,fval,e] = fminconIPasQPsolver(H,f,A,b,mu_lower,mu_upper,mu_init)
%
% Solve a quadratic problem with the interior-point algorithm of fmincon
%load quadprog_error

opt = optimset( ...%'Display','iter', %
    'Hessian','user-supplied', ...
    'HessFcn',@(x,y) qphess(x,0,H), 'GradObj', 'on', ...
    'Algorithm','interior-point','TolX',10^-8);

[this_mu,fval,e] = fmincon(@(x) objective(x,H,f),mu_init,A,b,[],[],mu_lower,mu_upper,[],opt);

%[this_mu,fval,e] = quadprog(H, f, A, b, [], [], mu_lower, mu_upper, this_mu, opt);

%=================================================
function [fval,grad] = objective(x,H,f)
fval = 0.5*x'*H*x + f'*x;
grad = f;

%=================================================
function hess = qphess(x,lambda,H) %#ok x,lambda are not needed for QPs
hess = H;

