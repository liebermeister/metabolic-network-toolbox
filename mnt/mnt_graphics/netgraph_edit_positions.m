function network = netgraph_edit_positions(network,table_positions,flag_KEGG_ids,gp,reaction_KEGG_IDs)

% network = netgraph_edit_positions(network,table_positions,flag_KEGG_ids)
% flag_KEGG_ids  which name field in 'network' to use: metabolite_KEGGID or metabolites

eval(default('flag_KEGG_ids','0','gp','struct','reaction_KEGG_IDs','[]'));

if length(reaction_KEGG_IDs),
  network.reaction_KEGGID = reaction_KEGG_IDs;
end

gp_default.actprintnames = 0;
gp = join_struct(gp_default,gp);

network = netgraph_read_positions(network,table_positions,[0,0],1,flag_KEGG_ids);
network = netgraph_move(network,'all',gp);
netgraph_print_positions(network,table_positions,[0,0],'replace elements',flag_KEGG_ids);

display(sprintf('New positions saved to file %s',table_positions));