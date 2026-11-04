function netgraph_wheelplot(network,network_CoHid,data,options,wheelpars)

% NETGRAPH_WHEELPLOT Draw a network with a wheel chart on each reaction
%
% netgraph_wheelplot(network,network_CoHid,data,options,wheelpars)
%
% network        network whose reaction positions are used for the wheel charts
% network_CoHid  network that is drawn (e.g. without cofactors)
% data           (#reactions x #segments matrix) values for the wheel charts
% options        graphics options for network_CoHid (struct; see 'help netgraph_draw')
% wheelpars      options for the wheel charts (see wheelchart); default radius 0.03

[nm,nr] = size(network.N);

netgraph_concentrations(network_CoHid,[],[],1,options);

xy = network.graphics_par.x(:,nm+(1:nr));

if ~isfield(wheelpars,'radius'), wheelpars.radius = 0.03; end

for it = 1:nr
  if isfinite(data(it,1))
    wheelchart(data(it,:),wheelpars,xy(1,it),xy(2,it),wheelpars.radius);
  end
end
