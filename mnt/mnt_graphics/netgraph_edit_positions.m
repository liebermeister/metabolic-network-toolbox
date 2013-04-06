function network = netgraph_edit_positions(network,table_positions,flag_KEGG_ids)

% network = netgraph_edit_positions(network,table_positions,flag_KEGG_ids)

eval(default('flag_KEGG_ids','0'));

network = netgraph_read_positions(network,table_positions,[0,0],1,flag_KEGG_ids);
network = netgraph_move(network,'all',struct('actprintnames',0));
netgraph_print_positions(network,table_positions,[0,0],'replace elements',flag_KEGG_ids);

display(sprintf('New positions save to file %s',table_positions));