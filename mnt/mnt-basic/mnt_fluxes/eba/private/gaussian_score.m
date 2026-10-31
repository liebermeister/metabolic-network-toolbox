function [score,gradient] = gaussian_score(x,x_mean,x_cov,x_cov_inv)

% score = gaussian_score(x,x_mean,x_cov)
%
% Score function of a gaussian function (-log likelihood, up to a constant)
% for multivariate gaussian distribution
% x can be a matrix (many column vectors to be scored)

eval(default('x_cov_inv','pinv(x_cov)'));

for it = 1:size(x,2),
  score(it)      = 0.5 * (x(:,it)-x_mean)' * x_cov_inv * (x(:,it)-x_mean);
end

if nargout>1,
  for it = 1:size(x,2),
    gradient(:,it) = x_cov_inv * (x(:,it)-x_mean);
  end
end
