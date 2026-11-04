function netgraph_fluxes(network,J,varargin)

% NETGRAPH_FLUXES Draw a network with fluxes
%
% netgraph_fluxes(network, J, options)
% netgraph_fluxes(network, J, 'name1', value1, ...)
%
% J   fluxes: shown by arrows, and by reaction symbols (absolute values)
% options graphics options (struct or name-value pairs; see 'help netgraph_draw')
%
% See also: netgraph_concentrations, netgraph_draw

options = netgraph_options(varargin, network);

options_default = struct('actvalues',abs(J),'arrowvalues',J,'arrowvaluesmax',max(abs(J)),'arrowstyle','fluxes');

options = join_struct(options_default,options);

netgraph_draw(network,options);
