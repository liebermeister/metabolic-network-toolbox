function synergy_network_plot(network, reaction_values, reaction_pair_values, col, varargin)

% synergy_network_plot(network, reaction_values, reaction_pair_values, col, gp)
%
% Arguments:
% reaction_values       (vector) values for single reactions
% reaction_pair_values  (matrix) values for pairs of reactions
%                                  only upper triangle values are shown!
% col                   colormap
% gp                    graphics options (struct or name-value pairs; see 'help netgraph_draw'), plus:
%                          gp.show_diagonal_values
%                          gp.normalise_values
%                          gp.relative_threshold
%
% FOR METABOLITE SYNERGIES SEE netgraph_metabolite_synergies

eval(default('col','[]'));
gp = netgraph_options(varargin, network);

if isempty(col), col = rb_colors; end

gp_def = struct('actprintnames',0,'metprintnames',1,'arrowstyle','none','arrowsize',0.03, 'linecolor',[0 0 0],'arrowcolor',[.7 .7 .7],'colorbar',0, 'FontSize',8,'text_offset',[.01,-.01],'colormap',col,'hold_on',1,'linewidth',1, 'show_diagonal_values',0,'normalise_values',1,'flag_edges',1,'relative_threshold',0);

gp = join_struct(gp_def,gp);

[nm,nr] = size(network.N);

if isfield(network.graphics_par, 'reaction_mapping')
  actmap = network.graphics_par.reaction_mapping;
else
  actmap = 1:nr;
end

if gp.show_diagonal_values
  reaction_values = diag(reaction_pair_values);
  gp.actstyle = 'fixed';
end

% normalise by maximal absolute value
if gp.normalise_values
  if any(isfinite(reaction_values))
    reaction_values = reaction_values/max(abs(reaction_values(actmap)));
  end
  reaction_pair_values = reaction_pair_values - diag(diag(reaction_pair_values));
  reaction_pair_values = reaction_pair_values/max(max(abs(reaction_pair_values(actmap,actmap))));
end

if gp.relative_threshold
  reaction_pair_values(abs(reaction_pair_values)<gp.relative_threshold*max(abs(reaction_pair_values(:)))) = 0;
end

% the loop plots only upper triangle values!!
reaction_pair_values = reaction_pair_values';

% replace nan by zeros for graphics
reaction_pair_values(isnan(reaction_pair_values)) = 0;

[i1_list, i2_list, colors] = netgraph_pair_colors(reaction_pair_values, actmap, col);

if isfield(gp,'metvalues')
  c = gp.metvalues;
else
  c = [];
end

netgraph_concentrations(network,c,reaction_values,1,gp); hold on;

draw_reaction_arcs(network,colors,i1_list,i2_list,gp.linewidth);  hold on;
