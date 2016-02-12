function [network, indices] = network_add_metabolites(network,new_metabolites,new_external);

% network = network_add_metabolites(network,new_metabolites,new_external);
%
% the new metabolites are added at teh end of the metabolite list
  
new_metabolites = column(setdiff(new_metabolites,network.metabolites));

[nr,nm] = network_numbers(network);
nn = length(new_metabolites);

indices = nm+[1:nn];

if ~exist('new_external','var'),
  new_external =  ones(nn,1);
end
network.metabolites       = [network.metabolites; new_metabolites];
network.external          = [network.external; new_external];
network.N                 = [network.N; zeros(nn,nr)];
network.regulation_matrix = [network.regulation_matrix, zeros(nr,nn)];

if isfield(network,'s_init'),
  network.s_init            = [network.s_init; nan*zeros(nn,1)];
end

if isfield(network,'metabolite_names'),
  network.metabolite_names = [network.metabolite_names; new_metabolites];
end

if isfield(network,'metabolite_compartments'),
  network.metabolite_compartments = [network.metabolite_compartments; repmat({''},nn,1)];
end

if isfield(network,'metabolite_KEGGID'),
  network.metabolite_KEGGID = [network.metabolite_KEGGID; repmat({'unknown_KEGG_ID'},nn,1)];
end
