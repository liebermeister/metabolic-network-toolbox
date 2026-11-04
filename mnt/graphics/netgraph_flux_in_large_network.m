% NETGRAPH_FLUX_IN_LARGE_NETWORK Show the reactions with the largest fluxes in a large network
%
% [nn, M, n_show] = netgraph_flux_in_large_network(network, v, perc_show, table_positions, cofactors, not_cofactors, options)
%
% For a large network and a flux distribution v, returns a smaller network nn for display:
% without cofactors, without metabolites with more than 5 reactions (sum of
% absolute stoichiometric coefficients), and without metabolites
% that take no part in the reactions with the largest fluxes (|v| >= perc_show * max |v|).
% To draw it with the fluxes, use netgraph_concentrations(nn, [], v, 1).
%
% Inputs
%   network          metabolic network
%   v                flux vector
%   perc_show        (optional, default 0.002) threshold for large fluxes, relative to the largest flux
%   table_positions  (optional) layout file or SBtab table; positions in nn are read from it
%   cofactors        (optional) list of metabolites to be omitted
%   not_cofactors    (optional) list of metabolites to be kept in any case
%   options          (optional) struct or name-value pairs, with field
%                      verbose  (default 1) plot the sorted absolute fluxes (in the current axes)
%                               and print the reactions with the largest fluxes
%
% Outputs
%   nn      network for display
%   M       (always empty; kept for compatibility)
%   n_show  number of reactions with large fluxes
%
% See also: netgraph_simple_graph, netgraph_concentrations

function [nn,M,n_show] = netgraph_flux_in_large_network(network,v,perc_show,table_positions,cofactors,not_cofactors,varargin)

eval(default('perc_show','0.002','table_positions','''''','cofactors','[]','not_cofactors','[]'));

options = join_struct(struct('verbose',1), netgraph_options(varargin, network));

v = v(:);
is_large = abs(v) >= perc_show * max(abs(v));
n_show   = sum(is_large);
M = [];

if options.verbose
  % distribution of flux values
  plot(sort(abs(v))); set(gca,'YScale','log'); xlabel('rank'); ylabel('absolute flux');
  title(sprintf('%d reactions above %g x maximum', n_show, perc_show));
  % reactions with the largest fluxes
  if isfield(network,'reaction_names'), names = network.reaction_names; else, names = network.actions; end
  [~,order] = sort(-abs(v));
  fprintf('Reactions with the %d largest fluxes\n', n_show);
  print_matrix(v(order(1:n_show)), names(order(1:n_show)));
end

% metabolites to be omitted: cofactors, metabolites in many reactions, and metabolites
% that take no part in the reactions with large fluxes
metabolite_rank = full(sum(abs(network.N), 2));   % sum of absolute stoichiometric coefficients
not_involved    = (abs(network.N) * double(is_large)) == 0;
omit = [cofactors(:); network.metabolites(metabolite_rank > 5); network.metabolites(not_involved)];
omit = setdiff(omit, not_cofactors);

nn = netgraph_simple_graph(network, omit);

if isfield(nn,'metabolite_names')
  nn.metabolites = nn.metabolite_names;
end
if ~isempty(table_positions)
  nn = netgraph_read_positions(nn, table_positions, [0,0], 1, 0);
end
if isfield(nn,'metabolite_names') && isfield(nn,'graphics_par')
  nn.graphics_par.metnames = nn.metabolite_names;
end
