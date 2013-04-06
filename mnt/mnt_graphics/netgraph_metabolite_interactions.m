% netgraph_metabolite_interactions(network, influence_values, interaction_values, col, gp)
%
% For reaction interactions, see 'interaction_network_plot'

function netgraph_metabolite_interactions(network, influence_values, interaction_values, col, gp)


% to make this script nicer, update it according to interaction_network_plot

eval(default('col','my_colors(250)','gp','struct'));

plot_parameters = struct('actprintnames',0,'metprintnames',1,'arrowsize',0.03,...
    'linecolor',[0 0 0],'arrowcolor',[.7 .7 .7],'colorbar',0,...
    'FontSize',8,'text_offset',[.01,-.01],'colormap',col,'hold_on',1);

plot_parameters = join_struct(plot_parameters,gp);

[nm,nr] = size(network.N);

metmap = network.graphics_par.metabolite_mapping;

if size(interaction_values,1) ~= max(metmap), error('Wrong matrix size'); end

for i1 = 1:nm,
  for i2 = 1:i1,
    colors(:,i1,i2) = [nan;nan;nan];
    if isfinite(interaction_values(metmap(i1),metmap(i2))),
      if sign(interaction_values(metmap(i1),metmap(i2))),
        colors(:,i1,i2) = col(1+floor((size(col,1)-1)*(interaction_values(metmap(i1),metmap(i2))+1)/2),:);
      end
    end
  end;
end

netgraph_concentrations(network,influence_values,[],1,plot_parameters); hold on;
draw_metabolite_arcs(network,colors,2);   
netgraph_concentrations(network,influence_values,[],1,plot_parameters); hold off;
