function network = netgraph_read_positions(network,table_positions,offsets,fill_nans,flag_KEGG_ids,reaction_KEGG_IDs)

% network = netgraph_read_positions(network,table_positions,offsets,fill_nans,flag_KEGG_ids)
%
% flag_KEGG_ids  which name field in 'network' to use: metabolite_KEGGID or metabolites

eval(default('offsets','[]','fill_nans','1','flag_KEGG_ids','0','reaction_KEGG_IDs','[]'));

if isempty(offsets),
  offsets = [0,0];
end

x     = nan*ones(2,length(network.metabolites)+length(network.actions));
fixed = zeros(1,length(network.metabolites)+length(network.actions));

if exist('sbtab_version','file'), 
  sbtab_table = sbtab_table_load(table_positions);
  columns = sbtab_table_get_all_columns(sbtab_table);
  A{1} = columns.Element;
  A{2} = cell_string2num(columns.PositionX);
  A{3} = cell_string2num(columns.PositionY);
else,
  fid = fopen(table_positions);
  A   = textscan(fid,'%s%f%f','delimiter','\t');
  fclose(fid);
end

% ---

switch flag_KEGG_ids,
  case 0,
    model_metabolite_names = network.metabolites;
    model_reaction_names   = network.actions;
  case 1,
    model_metabolite_names = network.metabolite_KEGGID;
    if length(reaction_KEGG_IDs),
      model_reaction_names = reaction_KEGG_IDs;
    else
      model_reaction_names   = network.reaction_KEGGID;
    end
end

if length(unique(model_metabolite_names)) <  length(model_metabolite_names),
  error('non-unique metabolite names');
end

if length(unique(model_reaction_names)) <  length(model_reaction_names),
  error('non-unique reaction names');
end

l1 = label_names(A{1},model_metabolite_names);
l1(isnan(A{2})) = 0;
ok = find(l1);

x(1,l1(ok)) = A{2}(ok)';
x(2,l1(ok)) = A{3}(ok)';
fixed(find(label_names(model_metabolite_names,A{1}))) = 1;

% ---

l1 = label_names(A{1},model_reaction_names);
l1(isnan(A{2})) = 0;
ok = find(l1);

x(1,length(network.metabolites)+l1(ok)) = A{2}(ok)';
x(2,length(network.metabolites)+l1(ok)) = A{3}(ok)';
fixed(length(network.metabolites)+find(label_names(model_reaction_names,A{1}))) = 1;

if fill_nans,
  x(:,find(sum(isnan(x)))) = repmat([nanmean(x(1,:)'); nanmean(x(2,:)')],1,length(find(sum(isnan(x)))));
  x(find(isnan(x))) = randn(size(x(find(isnan(x)))));
end

x = x - repmat(offsets',1,size(x,2));

if ~isfield(network,'graphics_par'),
  network = netgraph_make_graph(network,[],[]);
end
network.graphics_par.x = x;
network.graphics_par.fixed = fixed;
