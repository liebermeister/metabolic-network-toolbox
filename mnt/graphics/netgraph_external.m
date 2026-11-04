% NETGRAPH_EXTERNAL Mark external metabolites with the mouse
%
% [network, p] = netgraph_external(network, draw_only, external)
%
% Mouse buttons: left/centre: mark/unmark a metabolite; right: finish
%
% Inputs:
%   network    network structure; network.external contains the current choice
%   draw_only  (optional, default 0) only show the current choice
%   external   (optional) initial choice (replaces network.external)
%
% Outputs:
%   network  network with the new field 'external'
%   p        graphics settings used
%
% See also: netgraph_mark_metabolites

function [network,p] = netgraph_external(network,draw_only,external)

if exist('external','var'),   network.external = external; end
if ~exist('draw_only','var'), draw_only = 0; end

[network.external,p] = netgraph_mark_metabolites(network,network.external,draw_only);

title('External metabolites confirmed');
