% NETWORK_FIND_METABOLITES Find the indices of metabolites by their names
%
%indices = network_find_metabolites(network,metabolite)
%
%INPUT
%metabolites: metabolite name or list of metabolite names
%
% See also: network_find_reactions

function ind = network_find_metabolites(network,metabolite)

if ischar(metabolite), metabolite={metabolite}; end

ind = label_names(metabolite,network.metabolites,'single');
