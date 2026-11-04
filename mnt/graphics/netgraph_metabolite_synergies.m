% NETGRAPH_METABOLITE_SYNERGIES Draw values for pairs of metabolites as arcs
%
% netgraph_metabolite_synergies(network, influence_values, synergy_values, col, gp)
% netgraph_metabolite_synergies(network, influence_values, synergy_values, col, 'name1', value1, ...)
%
% influence_values  (vector) values for single metabolites
% synergy_values    (matrix) values for pairs of metabolites
% col               colormap (default my_colors(250))
% gp                graphics options (struct or name-value pairs; see 'help netgraph_draw'), plus
%                     relative_threshold, normalise_values, show_diagonal_values (see synergy_network_plot)
%
% For reaction synergies, see 'synergy_network_plot'

function netgraph_metabolite_synergies(network, influence_values, synergy_values, col, varargin)

eval(default('col','my_colors(250)','influence_values','[]'));
gp = netgraph_options(varargin, network);

gp = join_struct(struct('relative_threshold',0,'normalise_values',1,'show_diagonal_values',0,'linewidth',1),gp);

plot_parameters = struct('actprintnames',0,'metprintnames',1,'arrowsize',0.03,...
    'linecolor',[0 0 0],'arrowcolor',[.7 .7 .7],'colorbar',0,...
    'FontSize',8,'text_offset',[.01,-.01],'colormap',col,'hold_on',1);
plot_parameters = join_struct(plot_parameters,gp);

[nm,nr] = size(network.N);

if isfield(network.graphics_par, 'metabolite_mapping')
  metmap = network.graphics_par.metabolite_mapping;
else
  metmap = 1:nm;
end

if size(synergy_values,1) < max(metmap)
  error('Wrong matrix size');
end

if gp.show_diagonal_values
  influence_values = diag(synergy_values);
end

% normalise by maximal absolute value
if gp.normalise_values
  if any(isfinite(influence_values))
    influence_values = influence_values/max(abs(influence_values(metmap)));
  end
  synergy_values = synergy_values - diag(diag(synergy_values));
  synergy_values = synergy_values/max(max(abs(synergy_values(metmap,metmap))));
end

if gp.relative_threshold
  synergy_values(abs(synergy_values)<gp.relative_threshold*max(abs(synergy_values(:)))) = 0;
end

% the loop plots only upper triangle values!!
synergy_values = synergy_values';

% replace nan by zeros for graphics
synergy_values(isnan(synergy_values)) = 0;

[i1_list, i2_list, colors] = netgraph_pair_colors(synergy_values, metmap, col);

netgraph_concentrations(network,influence_values,[],1,plot_parameters); hold on;

draw_metabolite_arcs(network,colors,i1_list,i2_list,gp.linewidth);
