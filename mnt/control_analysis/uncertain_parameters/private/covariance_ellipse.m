function ellipse = covariance_ellipse(data)

m     = mean(data,2);
[v,d] = eig( cov(data') );

x=-1:0.01:1;
y1=sqrt(1-x.^2);
y2=-y1;

circle = [[x,fliplr(x)];[y1,y2]];
%plot(data(1,:),data(2,:),'.'); hold on;

ellipse = v' * sqrtm(d) * v * circle + diag(m)*ones(size(circle));

%plot(ellipse(1,:),ellipse(2,:));
