%indices = network_find_metabolites(network,metabolites)
%
%INPUT
%metabolites: metabolite name or list of metabolite names

function ind = network_find_metabolites(network,metabolite)

if isstr(metabolite), metabolite={metabolite}; end
ind=label_names(metabolite,network.metabolites,'single');
