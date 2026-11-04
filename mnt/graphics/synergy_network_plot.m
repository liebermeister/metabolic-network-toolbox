function synergy_network_plot(network, reaction_values, reaction_pair_values, col, varargin)

% SYNERGY_NETWORK_PLOT Draw values for pairs of reactions as arcs between the reactions
%
% synergy_network_plot(network, reaction_values, reaction_pair_values, col, options)
%
% Inputs:
%   network               network structure with layout (see netgraph_make_graph)
%   reaction_values       (vector) values for single reactions (shown by the reaction symbols)
%   reaction_pair_values  (matrix) values for pairs of reactions (only the upper triangle is shown)
%   col                   (optional, default rb_colors) colormap
%   options               graphics options (struct or name-value pairs; see 'help netgraph_draw'), plus
%                           show_diagonal_values  (default 0) show the diagonal of reaction_pair_values
%                                                 as reaction values
%                           normalise_values      (default 1) scale values by their maximal absolute value
%                           relative_threshold    (default 0) omit pair values below this fraction of the maximum
%
% See also: netgraph_metabolite_synergies, netgraph_draw

eval(default('col','[]'));
options = netgraph_options(varargin, network);

if isempty(col), col = rb_colors; end

options_def = struct('actprintnames',0,'metprintnames',1,'arrowstyle','none','arrowsize',0.03, 'linecolor',[0 0 0],'arrowcolor',[.7 .7 .7],'colorbar',0, 'FontSize',8,'text_offset',[.01,-.01],'colormap',col,'hold_on',1,'linewidth',1, 'show_diagonal_values',0,'normalise_values',1,'flag_edges',1,'relative_threshold',0);

options = join_struct(options_def,options);

[nm,nr] = size(network.N);

if isfield(network.graphics_par, 'reaction_mapping')
  actmap = network.graphics_par.reaction_mapping;
else
  actmap = 1:nr;
end

if options.show_diagonal_values
  reaction_values = diag(reaction_pair_values);
  options.actstyle = 'fixed';
end

% normalise by maximal absolute value
if options.normalise_values
  if any(isfinite(reaction_values))
    reaction_values = reaction_values/max(abs(reaction_values(actmap)));
  end
  reaction_pair_values = reaction_pair_values - diag(diag(reaction_pair_values));
  reaction_pair_values = reaction_pair_values/max(max(abs(reaction_pair_values(actmap,actmap))));
end

if options.relative_threshold
  reaction_pair_values(abs(reaction_pair_values)<options.relative_threshold*max(abs(reaction_pair_values(:)))) = 0;
end

% the loop plots only upper triangle values!!
reaction_pair_values = reaction_pair_values';

% replace nan by zeros for graphics
reaction_pair_values(isnan(reaction_pair_values)) = 0;

[i1_list, i2_list, colors] = netgraph_pair_colors(reaction_pair_values, actmap, col);

if isfield(options,'metvalues')
  c = options.metvalues;
else
  c = [];
end

netgraph_concentrations(network,c,reaction_values,1,options); hold on;

draw_reaction_arcs(network,colors,i1_list,i2_list,options.linewidth);  hold on;
