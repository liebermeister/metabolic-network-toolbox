function networkG = netgraph_layout_from_table(network, layout_table)

% NETGRAPH_LAYOUT_FROM_TABLE Network for graphics from a network and an SBtab layout table
%
% networkG = netgraph_layout_from_table(network, layout_table)
%
% Rebuild the network for graphics (networkG) from a model network and an SBtab layout table
% (TableID 'NetworkLayout', written by netgraph_layout_table): metabolites and reactions that are not
% listed in the table are removed (netgraph_simple_graph), and the positions and invisible flags
% are read from the table (netgraph_read_positions).
%
% (Until 2026-11-17 the IMF function imf_layout_network_from_table.)
%
% See also: netgraph_layout_table, mnt_load_model_and_state

elements = sbtab_table_get_column(layout_table, 'Element');

hide_metabolites = network.metabolites(~ismember(network.metabolites, elements));
hide_reactions   = find(~ismember(network.actions, elements));

networkG = netgraph_make_graph(network);
networkG = netgraph_simple_graph(networkG, hide_metabolites, [], [], hide_reactions);
networkG = netgraph_read_positions(networkG, layout_table, [], 0);
