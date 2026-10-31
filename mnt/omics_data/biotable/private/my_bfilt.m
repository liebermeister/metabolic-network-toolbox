function r = my_bfilt(data)

% r = my_bfilt(data)
%
% Smooth data along rows (weighted moving average over 5 neighbouring values,
% omitting the largest and smallest of them); the data need at least 5 columns

clear r;

for it =1:size(data,1)

  dum2 = [];

  for it2 =3:size(data,2)-2
    dum = data(it,it2-2:it2+2);
    dum(find(~isfinite(dum)))= mean(dum','omitnan');
    [d1,order] = sort(dum);
    if d1(1)==d1(end)
      dum2 = [dum2; d1(1:3)];
    else
      dum([order(1),order(end)])  = nan;
      dum2 = [dum2; dum(isfinite(dum))];
    end
  end
  r(it,:) = (dum2*[0.25 0.5 0.25]')';
end

r = [data(:,1:2)*[0.7,0.3]',  data(:,1:3)*[0.25 0.5 0.25]', r ,...
     data(:,end-2:end)*[0.25 0.5 0.25]', data(:,end-1:end)*[0.3,0.7]'];
