% output = netgraph_print_positions(network,filename,offsets,policy,flag_KEGG_ids)

function output = netgraph_print_positions(network,filename,offsets,policy,flag_KEGG_ids)

eval(default('filename','[]','offsets','[0,0]','policy','''replace elements''','flag_KEGG_ids','0'));

[nr,nm] = network_numbers(network);

switch flag_KEGG_ids,
  case 0,
    metnames = network.metabolites;
    actnames   = network.actions;
  case 1,
    metnames = network.metabolite_KEGGID;
    actnames   = network.reaction_KEGGID;
end
metpositions = network.graphics_par.x(:,1:nm) + repmat(offsets',1,nm);

actpositions = network.graphics_par.x(:,nm+1:end) + repmat(offsets',1,nr);

names     = [metnames; actnames];
positions = [metpositions, actpositions];

switch policy, 
  case 'add nonexistent',  
    t         = sbtab_table_load(filename);
    elements  = sbtab_table_get_column(t,'Element');
    x         = cell_string2num(sbtab_table_get_column(t,'PositionX'));
    y         = cell_string2num(sbtab_table_get_column(t,'PositionY'));
    ll        = label_names(names,elements); ind_keep = find(ll==0);
    names     = [elements; names(ind_keep)];
    positions = [[x';y'], positions(:,ind_keep) ];
  case 'replace elements',  
    t = sbtab_table_load(filename);
    elements  = sbtab_table_get_column(t,'Element');
    x         = cell_string2num(sbtab_table_get_column(t,'PositionX'));
    y         = cell_string2num(sbtab_table_get_column(t,'PositionY'));
    ll        = label_names(elements,names); ind_keep = find(ll==0);
    names     = [elements(ind_keep); names];
    positions = [ [x(ind_keep)';y(ind_keep)'], positions];
  case 'replace file',  
end

fid = fopen(filename,'w');

fprintf(fid,'!Element\t!PositionX\t!PositionY\n');

for it = 1:length(names),
  fprintf(fid,'%s\t%f\t%f\n',names{it},positions(1,it),positions(2,it));
end
fclose(fid);
%output = table([names,cellstr(num2str(positions(1,:)')),cellstr(num2str(positions(2,:)'))],0,filename);