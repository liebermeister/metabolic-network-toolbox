function draw_reaction_arcs(network, colors, linewidth)

% draw_reaction_arcs(network,colors,linewidth)

eval(default('linewidth','2'));

[nm,nr] = size(network.N);

xcentre = mean(network.graphics_par.x(:,nm+1:end),2);

hold on
for i1 = 1:nr,
  for i2 = 2:i1,
   x1 = network.graphics_par.x(:,nm+i1);
   x2 = network.graphics_par.x(:,nm+i2);
   arcsign = det([((x1+x2)/2 -xcentre), x1-x2]) >=0;
   x = arc(x1,x2, 2*(arcsign-0.5)*.2);
   col = squeeze(colors(:,i1,i2));
   if sum(isnan(col))==0,
     plot(x(1,:),x(2,:),'color',col,'Linewidth',linewidth); 
   end
  end
end
