% Graphical display and manipulation of metabolic networks
%
% Layout and display settings of a network are stored in the field 'graphics_par'
% (created by netgraph_make_graph; see 'help metabolic_network').
% For all display options, see 'help netgraph_draw'.
%
% Example:
%   network = netgraph_make_graph(network);              % compute a layout
%   netgraph_concentrations(network, c, v, 1);           % show concentrations c and fluxes v
%   network = netgraph_move(network);                    % move nodes with the mouse
%   netgraph_print_positions(network, 'layout.tsv');     % save the layout
%   network = netgraph_read_positions(network, 'layout.tsv');  % read it again later
%
% ---------------------------------------------------------------------------
% Draw networks
%   netgraph_draw                  draw a network (all display options are listed here)
%   netgraph_concentrations        draw a network with metabolite levels and fluxes
%   netgraph_fluxes                draw a network with fluxes
%   netgraph_regulation            draw a network with regulation arrows (from network.kinetics)
%   netgraph_draw_simple           draw nodes and lines only
%   netgraph_draw_connections      draw lines between reactions (from a matrix of relations)
%   netgraph_wheelplot             draw a network with a wheel chart on each reaction
%   synergy_network_plot           draw values for pairs of reactions as arcs
%   netgraph_metabolite_synergies  draw values for pairs of metabolites as arcs
%   draw_reaction_arcs             draw arcs between reactions
%   draw_metabolite_arcs           draw arcs between metabolites
%   netgraph_mark_reaction         draw a circle around a reaction
%   netgraph_flux_in_large_network show the reactions with the largest fluxes in a large network
%
% Movies
%   netgraph_movie                 movie of metabolite levels or fluxes over time
%   netgraph_double_movie          movie of two data sets at the same time
%   netgraph_flux_movie            movie of one flux distribution (moving arrows)
%   netgraph_flux_arrows_movie     movie of several flux distributions
%
% Layout
%   netgraph_make_graph            create the field 'graphics_par' and compute a layout
%   netgraph_calculate_positions   compute node positions (force-directed layout)
%   netgraph_find_x                compute node positions (used by netgraph_calculate_positions)
%   netgraph_move                  move nodes with the mouse
%   netgraph_edit_positions        read a layout file, move nodes with the mouse, and save the file
%   netgraph_read_positions        read a layout file (SBtab table)
%   netgraph_print_positions       save a layout file (SBtab table)
%   netgraph_assign_positions      choose node positions on a picture with the mouse
%   netgraph_save_positions        save positions in an older file format (deprecated)
%   netgraph_clear                 reset display settings, keeping the layout
%
% Networks for display
%   netgraph_simple_graph          network without some metabolites (e.g. cofactors)
%   netgraph_split_metabolites     network in which some metabolites appear once per reaction
%   netgraph_cofactor_graphs       both of the above, for cofactors
%
% Choose elements with the mouse
%   netgraph_mark_metabolites      mark metabolites
%   netgraph_mark_reactions        mark reactions
%   netgraph_external              mark external metabolites
%   netgraph_reversible            mark reversible reactions
%   netgraph_define_module         mark metabolites that form a module
%   netgraph_remove                remove metabolites and reactions
