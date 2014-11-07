function [network_CoHid,network_CoSplit] = netgraph_cofactor_graphs(network, cofactors, table_positions, max_met_degree,no_cofactors)

% [network_CoHid,network_CoSplit] = netgraph_cofactor_graphs(network, cofactors, table_positions, max_met_degree)

eval(default('cofactors','[]','table_positions','[]','max_met_degree','5','no_cofactors','[]'));

if length(max_met_degree),
  met_degree = sum(abs(sign(network.N)),2);
  cofactors  = unique([column(cofactors); network.metabolites(met_degree>max_met_degree)]);
end

cofactors  = setdiff(cofactors,no_cofactors);

network_CoHid   = netgraph_simple_graph(network,cofactors);
if nargout >1,
  network_CoSplit = netgraph_split_metabolites(network,cofactors);
end

if length(table_positions),
  network_CoHid   = netgraph_read_positions(network_CoHid,table_positions,[0,0],1);
  if nargout >1,
    network_CoSplit = netgraph_read_positions(network_CoSplit,table_positions,[0,0],1);
  end
end

%display('Please wait ...');
%pars = struct('initial_relax',200,'manual',0,'actprintnames',0,'arrowstyle','directions','arrowsize',0.03);
%network_CoSplit = netgraph_move(network_CoSplit,'text',pars);
%network_CoHid   = netgraph_move(network_CoSplit,'text',pars);
