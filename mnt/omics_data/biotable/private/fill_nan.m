function M = fill_nan(M)

% M = fill_nan(M)
%
% In each row of M, replace missing values (NaN) before the first and after
% the last finite value by these values

for it = 1:size(M,1)
  dum = find(isfinite(M(it,:)));
  if isempty(dum), continue; end
  M(it,1:min(dum)-1)   = M(it,min(dum));
  M(it,max(dum)+1:end) = M(it,max(dum));
end
