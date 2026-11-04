% NETGRAPH_REVERSIBLE Mark reversible reactions with the mouse
%
% [network, p] = netgraph_reversible(network)
%
% Mouse buttons: left/centre: mark/unmark a reaction; right: finish.
% The choice is stored in network.reversible.
%
% See also: netgraph_mark_reactions

function [network,p] = netgraph_reversible(network)

[network.reversible,p] = netgraph_mark_reactions(network,network.reversible);

title('Reversible reactions confirmed');
