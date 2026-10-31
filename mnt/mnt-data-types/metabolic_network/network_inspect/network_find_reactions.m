% NETWORK_FIND_REACTIONS Find the reactions in which metabolites participate
%
%indices = network_find_reactions(network,metabolites)
%
%find reactions in which a certain metabolite participates
%metabolites: metabolite name or list of metabolite names
%
% See also: network_find_metabolites

function re = network_find_reactions(network,metabolite)

only_one=0;
if ischar(metabolite), metabolite={metabolite}; only_one=1; end
ind = label_names(metabolite,network.metabolites,'single');

if only_one
  re=find(network.N(ind,:)~=0);
else
  for i=1:length(ind)
    re{i}=find(network.N(ind(i),:)~=0);
  end
end
