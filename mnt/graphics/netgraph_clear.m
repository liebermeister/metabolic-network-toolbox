% NETGRAPH_CLEAR Reset the display settings of a network, keeping the layout
%
% network = netgraph_clear(network)
%
% Replaces network.graphics_par by default display settings; node positions (x) and
% the adjacency matrices (m, db) are kept.
%
% See also: netgraph_make_graph, netgraph_draw

function network = netgraph_clear(network)

p.metnames = network.metabolites;
p.actnames = network.actions;

p.metstyle       = 'box';
p.metcolstyle    = 'values';
p.metcol         = 'c';
p.metprintnames  = 1;
p.metprintvalues = 0;
p.metvalues      = zeros(length(network.metabolites),1);
p.metvalues_std  = zeros(length(network.metabolites),1);

p.actstyle       = 'box';
p.actcolstyle    = 'values';
p.actcol         = 'r';
p.actprintnames  = 1;
p.actprintvalues = 0;
p.actvalues      = zeros(length(network.actions),1);
p.actvalues_std  = zeros(length(network.actions),1);

p.arrowstyle     = 'none';
p.arrowvalues    = p.actvalues;

p.squaresize = 1;
p.arrowsize  = 0.5;
p.lambda     = 1;

p.x          = network.graphics_par.x;
p.m          = network.graphics_par.m;
p.db          = network.graphics_par.db;

network.graphics_par = p;
