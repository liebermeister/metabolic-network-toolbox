function [mu,sig]=resampling(data,n);

% calculate bootstrap samples of mean and standard deviation [mu,sig]=resampling(data,n);

mu=[];
sig=[];
for it=1:n,
  ind = ceil(length(data)*rand(length(data),1));
  mem = data(ind);
  mu=[mu mean(mem)];
  sig=[sig std(mem)];
end
