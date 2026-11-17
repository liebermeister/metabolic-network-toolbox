function networkG = netgraph_make_layout_network(network, hide_metabolites, layout_file, hide_reactions)

% NETGRAPH_MAKE_LAYOUT_NETWORK Network for graphics, with a layout from a file and hidden elements
%
% networkG = netgraph_make_layout_network(network, hide_metabolites, layout_file, hide_reactions)
%
% Make an extra network struct for graphics (e.g. for netgraph_concentrations), in which some
% metabolites (e.g. cofactors) and reactions are not shown, with positions from a layout file.
%
%   network            network struct
%   hide_metabolites   metabolite IDs (network.metabolites) not shown, e.g. cofactors
%   layout_file        SBtab file (or SBtab table) with the positions of metabolites and reactions
%                      (table Position; see netgraph_read_positions)
%   hide_reactions     (optional) reaction IDs (network.actions) not shown (default: the biomass
%                      reactions network.biomass_reactions, if this field exists; otherwise {})
%
% The layout of networkG can be stored in an SBtab model file with netgraph_layout_table and
% rebuilt with netgraph_layout_from_table (see mnt_save_model_and_state, mnt_load_model_and_state).
%
% (Until 2026-11-17 the IMF function imf_make_layout_network.)
%
% See also: netgraph_simple_graph, netgraph_read_positions, netgraph_layout_table

eval(default('hide_metabolites','{}','hide_reactions','[]'));
if isempty(hide_reactions), hide_reactions = network_biomass_reactions(network); end
if ischar(hide_reactions), hide_reactions = {hide_reactions}; end

networkG = netgraph_make_graph(network);
networkG = netgraph_simple_graph(networkG, hide_metabolites);
networkG = netgraph_read_positions(networkG, layout_file);

ind_hide = label_names(hide_reactions, networkG.actions);
if any(ind_hide == 0)
  warning('netgraph_make_layout_network: reaction(s) %s not found in the network', strjoin(hide_reactions(ind_hide == 0), ', '));
end
networkG.graphics_par.actinvisible(ind_hide(ind_hide > 0)) = 1;
