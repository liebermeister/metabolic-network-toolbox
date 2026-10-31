function [m_mean,s_mean,CV_mean,m_std,s_std,CV_std] = simple_bootstrap(data,nit)

% [m_mean,s_mean,m_std,s_std] = simple_bootstrap(data,nit)
%
% data m x n matrix, n samples (m variables)
% nit : number of boostrap samples
%
% m_mean  estimated mean value
% s_mean  estimated standard deviation
% m_std   error of mean value
% s_std   error of standard deviation

mem_mean = [];
mem_std  = [];
mem_CV   = [];

if ~exist('nit','var'), nit=100; end

for i=1:nit,
  sample   = data(:,ceil(size(data,2) * rand(1,size(data,2))));
  mem_mean = [mem_mean, mean(sample,2) ];
  mem_std  = [mem_std,  std( sample' )' ];
  mem_CV   = [mem_CV,  std( sample' )'./mean(sample,2) ];
end

m_mean  = mean(mem_mean,2);
s_mean  = mean(mem_std,2);
CV_mean = mean(mem_CV,2);
m_std   = std(mem_mean')';
s_std   = std(mem_std')';
CV_std  = std(mem_CV')';
