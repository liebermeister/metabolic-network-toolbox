function interaction_network_plot(network, reaction_values, reaction_pair_values, col, gp)

% interaction_network_plot(network, reaction_values, reaction_pair_values, col, gp)
%
% Arguments:
% reaction_values       (vector) values for single reactions 
% reaction_pair_values  (matrix) values for pairs of reactions 
%                                  only upper triangle values are shown!
% col                   colormap
% gp                    graphics parameters for network, but with extra entries: 
%                          gp.show_diagonal_values
%                          gp.normalise_values
%                          gp.relative_threshold
% 
% FOR METABOLITE INTERACTIONS SEE netgraph_metabolite_interactions

eval(default('col','[]','gp','struct'));

if isempty(col), col = rb_colors; end

gp_def = struct('actprintnames',0,'metprintnames',1,'arrowsize',0.03, 'linecolor',[0 0 0],'arrowcolor',[.7 .7 .7],'colorbar',0, 'FontSize',8,'text_offset',[.01,-.01],'colormap',col,'hold_on',1,'linewidth',2, 'show_diagonal_values',0,'normalise_values',1,'flag_edges',1,'relative_threshold',0);

gp = join_struct(gp_def,gp);

[nm,nr] = size(network.N);

if isfield(network.graphics_par, 'reaction_mapping'),
  actmap = network.graphics_par.reaction_mapping;
else,
  actmap = 1:nr;
end

%if length(actmap) ~= size(reaction_pair_values,1), error('Wrong matrix size'); end

if gp.show_diagonal_values,
  reaction_values = diag(reaction_pair_values);
end

if gp.normalise_values,
  if length(isfinite(reaction_values)),
    reaction_values = reaction_values/max(abs(reaction_values(actmap)));
  end
  reaction_pair_values = reaction_pair_values - diag(diag(reaction_pair_values));
  reaction_pair_values = reaction_pair_values/nanmax(nanmax(abs(reaction_pair_values(actmap,actmap))));
end

if gp.relative_threshold,
  reaction_pair_values(abs(reaction_pair_values)<gp.relative_threshold*max(abs(reaction_pair_values(:)))) = 0;
end

% the loop plots only upper triangle values!!
reaction_pair_values = reaction_pair_values';

% replace nan by zeros for graphics
reaction_pair_values(isnan(reaction_pair_values)) = 0;

for i1 = 1:nr, 
  for i2 = 2:i1,
    colors(:,i1,i2) = [nan;nan;nan];
    if isfinite(reaction_pair_values(actmap(i1),actmap(i2))),  
      if sign(reaction_pair_values(actmap(i1),actmap(i2))),  
        colors(:,i1,i2) = col(1+floor((size(col,1)-1)*(reaction_pair_values(actmap(i1),actmap(i2))+1)/2),:);
      end
    end
  end;
end

if isfield(gp,'metvalues'),
  c = gp.metvalues; 
else,
  c = [];
end

netgraph_concentrations(network,c,reaction_values,1,gp); hold on;

draw_reaction_arcs(network,colors,gp.linewidth);  hold on;
