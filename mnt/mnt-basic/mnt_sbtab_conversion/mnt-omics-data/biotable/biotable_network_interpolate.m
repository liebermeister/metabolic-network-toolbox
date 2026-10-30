function data = biotable_network_interpolate(data, network, type,recursive)

% data = biotable_network_interpolate(data, network, type, recursive)
%
% Replace missing values (NaN) in a biotable by the mean values of neighbouring
% elements in the network (type 'metabolites' or 'reactions')
% If recursive==1, repeat this as long as missing values can be filled

eval(default('recursive','0'));

data = biotable_clean(data);

[nm,nr] = size(network.N);
ne = size(data.DataMean,1);

NN = network.N;
NN(sum(abs(NN),2)>6,:) = 0; % disregard hub metabolites

switch type
  case 'metabolites'
    if nm ~= ne, error('wrong dimensions'); end
    M = double([double(NN~=0) * double(NN'~=0)]~=0);
  case 'reactions'
    if nr ~= ne, error('wrong dimensions'); end
    M = double([double(NN'~=0) * double(NN~=0)]~=0);
end

M = M-diag(diag(M));

for it = 1:ne
  ind = find(M(:,it));
  data_interpol(it,:) = mean(data.DataMean(ind,:),1,'omitnan');
end

n_missing = sum(isnan(data.DataMean(:)));
data.DataMean(isnan(data.DataMean)) =  data_interpol(isnan(data.DataMean));

% repeat as long as missing values can be filled
if recursive
  if sum(isnan(data.DataMean(:))) && sum(isnan(data.DataMean(:))) < n_missing
    data = biotable_network_interpolate(data, network, type,1);
  end
end
