function draw_metabolite_arcs(network, colors, linewidth)

% draw_metabolite_arcs(network, colors, linewidth)
%
%
% See also draw_reaction_arcs(network,colors,linewidth)

eval(default('linewidth','2'));

[nm,nr] = size(network.N);

xcentre = mean(network.graphics_par.x(:,1:nm),2);

hold on
for i1 = 1:nm,
  for i2 = 1:i1,
   x1 = network.graphics_par.x(:,i1);
   x2 = network.graphics_par.x(:,i2);
   arcsign = det([((x1+x2)/2 -xcentre), x1-x2]) >=0;
   x = arc(x1,x2, 2*(arcsign-0.5)*.2);
   col = squeeze(colors(:,i1,i2));
   if sum(isnan(col))==0,
   plot(x(1,:),x(2,:),'color',col,'Linewidth',linewidth); 
   end
  end
end
