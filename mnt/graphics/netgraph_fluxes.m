function netgraph_fluxes(network,J,varargin)

% NETGRAPH_FLUXES Draw a network with fluxes
%
% netgraph_fluxes(network, J, gp)
% netgraph_fluxes(network, J, 'name1', value1, ...)
%
% J   fluxes: shown by arrows, and by reaction symbols (absolute values)
% gp  graphics options (struct or name-value pairs; see 'help netgraph_draw')
%
% See also: netgraph_concentrations, netgraph_draw

gp = netgraph_options(varargin, network);

gp_default = struct('actvalues',abs(J),'arrowvalues',J,'arrowvaluesmax',max(abs(J)),'arrowstyle','fluxes');

gp = join_struct(gp_default,gp);

netgraph_draw(network,gp);
