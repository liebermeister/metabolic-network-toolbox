function new_network = network_add_reaction(network,reactants,stoich_coeffs,name,reversible,external); 

% new_network = network_add_reaction(network,reactants,stoich_coeffs,name,reversible); 

eval(default('reversible','1'));

reactants     = column(reactants);
stoich_coeffs = column(stoich_coeffs);

% remove metabolites that already exist in the network from the list "reactants"
ll = label_names(reactants,network.metabolites);
additional_metabolites = reactants(find(ll==0));

if length(additional_metabolites),
  additional_external = ones(size(additional_metabolites));
  nn = network_add_metabolites(network,additional_metabolites,additional_external);
else 
  nn = network;
end

[nr,nm] = network_numbers(nn);

new_network.metabolites               = nn.metabolites;
new_network.actions                   = nn.actions;
ind = label_names(reactants, new_network.metabolites);
new_network.actions                   = [network.actions; name];
new_network.external                  = nn.external;
new_network.N                         = nn.N;
new_network.N(:,nr+1)                 = 0;
new_network.N(ind,nr+1)               = stoich_coeffs;
new_network.regulation_matrix         = nn.regulation_matrix;
new_network.regulation_matrix(nr+1,:) = 0;
new_network.reversible                = [nn.reversible; reversible];

if isfield(network, 'metabolite_BIGGID'),
  new_network.metabolite_BIGGID = [network.metabolite_BIGGID; repmat({},length(additional_metabolites),1)];
end

if isfield(network, 'reaction_BIGGID'), new_network.reaction_BIGGID = [network.reaction_BIGGID; {[]}]; end
if isfield(network, 'enzyme_mass'),     new_network.enzyme_mass = [network.enzyme_mass; nan]; end

if isfield(network, 'reaction_names'),
  new_network.reaction_names = [network.reaction_names; {name}];
end

if isfield(network, 'formulae'),
  new_network.formulae = network_print_formulae(new_network);
end
