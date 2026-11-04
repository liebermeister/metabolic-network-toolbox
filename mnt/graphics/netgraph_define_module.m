% NETGRAPH_DEFINE_MODULE Choose the metabolites of a network module with the mouse
%
% [m, n] = netgraph_define_module(network)
%
% Mouse buttons: left/centre: mark/unmark a metabolite; right: finish
%
% Outputs:
%   m  indices of the chosen metabolites
%   n  (optional) subnetwork containing these metabolites (see network_choose)
%
% See also: netgraph_mark_metabolites

function [m,n] = netgraph_define_module(network)

[m,p] = netgraph_mark_metabolites(network);
m = find(m);

if nargout>1
  n = network_choose(network,m);
  n.kinetics = set_kinetics(n,network.kinetics.type);
end
