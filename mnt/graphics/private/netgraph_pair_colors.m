function [i1_list, i2_list, colors] = netgraph_pair_colors(pair_values, map, col)

% [i1_list, i2_list, colors] = netgraph_pair_colors(pair_values, map, col)
%
% Element pairs and colours for drawing pair values (e.g. synergies) as arcs
% (used by synergy_network_plot and netgraph_metabolite_synergies)
%
% pair_values  matrix of values (scaled to [-1,1]); only entries (map(i1),map(i2)) with i1 < i2 are used
% map          mapping from graph elements to rows/columns of pair_values
% col          colormap; values -1 and 1 are shown by the first and last colour
%
% Pairs with zero or non-finite values are omitted. Pairs are sorted by increasing absolute value,
% so the strongest arcs are drawn on top.

n = length(map);
i1_list    = [];
i2_list    = [];
value_list = [];

for i1 = 1:n
  for i2 = i1+1:n
    this_value = pair_values(map(i1),map(i2));
    if isfinite(this_value) && this_value ~= 0
      i1_list(end+1)    = i1;
      i2_list(end+1)    = i2;
      value_list(end+1) = this_value;
    end
  end
end

value_list = min(max(value_list,-1),1);
colors = col(1+floor((size(col,1)-1)*(value_list+1)/2),:);

[~,order] = sort(abs(value_list));
i1_list = i1_list(order);
i2_list = i2_list(order);
colors  = colors(order,:);
