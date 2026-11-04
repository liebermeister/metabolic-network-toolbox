% NETGRAPH_SIMPLE_GRAPH Network for display without some metabolites (e.g. cofactors)
%
% [n2, keep_met, keep_react] = netgraph_simple_graph(n1, names_met_hide, met_fixed, met_fixed_positions, ind_react_hide)
%
% Inputs:
%   n1                   network structure
%   names_met_hide       list of metabolite names to be removed
%   met_fixed            (optional) indices of metabolites with given positions
%   met_fixed_positions  (optional) positions of these metabolites (2 x n matrix)
%   ind_react_hide       (optional) indices of reactions to be removed
%
% Outputs:
%   n2          network without these metabolites and reactions, with a layout
%   keep_met    indices of the remaining metabolites (in n1)
%   keep_react  indices of the remaining reactions (in n1)
%
% See also: netgraph_split_metabolites, netgraph_cofactor_graphs

function [n2,keep_met,keep_react] = netgraph_simple_graph(n1,names_met_hide,met_fixed,met_fixed_positions,ind_react_hide)

eval(default('met_fixed','[]','met_fixed_positions','[]','ind_react_hide','[]'));

keep_met             = setdiff(1:length(n1.metabolites),network_find_metabolites(n1,names_met_hide));
%keep_react          = 1:length(n1.actions);
keep_react           = setdiff(find(sum(abs(n1.N(keep_met,:)),1)), ind_react_hide);
n2                   = network_subnetwork(n1,keep_met,keep_react);

if length(met_fixed)
  set_positions = network_find_metabolites(n2,met_fixed);
  x = zeros(2,length(n2.metabolites)+length(n2.actions));
  x(:,set_positions) = met_fixed_positions;
  n2 = netgraph_make_graph(n2,[],[],set_positions,x);
elseif isfield(n1,'graphics_par')
  [nm,nr] = size(n1.N);
  set_positions = network_find_metabolites(n1,n2.metabolites);
  x = [n1.graphics_par.x(:,set_positions), n1.graphics_par.x(:,nm+(1:nr))];
  n2 = netgraph_make_graph(n2,[],[],[1:length(n2.metabolites)+length(n2.actions)]',x);
else
  n2 = netgraph_make_graph(n2);
end

n2.graphics_par.metabolite_mapping = column(keep_met);
n2.graphics_par.reaction_mapping   = column(keep_react);
n2.graphics_par.metvalues = [];
n2.graphics_par.actvalues = [];

if isfield(n1,'graphics_par')
  if isfield(n1.graphics_par,'subplot_position')
    n2.graphics_par.subplot_position = n1.graphics_par.subplot_position;
  end
  if isfield(n1.graphics_par,'figure_position')
    n2.graphics_par.figure_position  = n1.graphics_par.figure_position;
  end
end
