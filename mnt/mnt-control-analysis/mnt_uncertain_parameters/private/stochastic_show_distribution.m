function stochastic_show_distribution(list,mu,sigma,positive,fignum)

% stochastic_show_distribution(list,mu,sigma,positive,fignum)
%
% list: monte carlo samples
% mu:  mean of logs (exact distribution)
% sigma: std dev  of logs (exact distribution)
% positive: (-1 or 1) for variables that are always negative
% fignum: figure number

if  exist('fignum','var'), figure(fignum);  end
if ~exist('positive','var'), positive=1;  end

list=list(find(abs(list)<10*abs(median(list))));

if length(list),
  nit = length(list);
  [n,bin] = hist(list,20);
  hist(list,20); hold on; colormap(1-gray)
  maxvalue = 1.1*max(bin);
  scale_y  = nit*(bin(2)-bin(1))
else,
  scale_y = 1;
%  [m1,s1] = lognormal_log2normal(mu,sigma);
  maxvalue = 3*exp(mu);
end

vals = (0:0.01:1)*maxvalue;
plot(vals*positive,scale_y*lognpdf(vals,mu,sigma),'LineWidth',2);
hold off;
