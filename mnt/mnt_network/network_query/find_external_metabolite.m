function index =  find_external_metabolite(network, metabolite)

% index =  find_external_metabolite(network, metabolite)

index = label_names({metabolite},network.metabolites(find(network.external)));