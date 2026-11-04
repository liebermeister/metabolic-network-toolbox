function network = netgraph_read_positions(network, layout_file, offsets, fill_nans, flag_KEGG_ids, reaction_KEGG_IDs, flag_fixed)

% NETGRAPH_READ_POSITIONS Read a network layout from a file
%
% network = netgraph_read_positions(network,layout_file,offsets,fill_nans,flag_KEGG_ids, reaction_KEGG_IDs, flag_fixed)
%
% Arguments:
%   'layout_file'   may either refer to a filename or to an SBtab data object
%   'fill_nans',    place elements missing in the layout file near their neighbours in the network
%                   (reactions between their metabolites, metabolites next to their reactions)
%   'flag_KEGG_ids' match elements (between model and layout file) by KEGG ids (default: 0)
%   'flag_fixed'    read information about fixed glyph positions from file
%
% If flag_KEGG_ids is set:
% use field 'network.metabolite_KEGGID' (or 'network.Identifiers_kegg_compound') for matching metabolites
% use argument list 'reaction_KEGG_IDs' or field 'network.reaction_KEGGID' for matching reactions
%
% See also: netgraph_print_positions, netgraph_edit_positions

eval(default('offsets','[]','fill_nans','1','flag_KEGG_ids','0','reaction_KEGG_IDs','[]','flag_fixed','1'));

if isempty(offsets), offsets = [0,0]; end

x = nan * ones(2,length(network.metabolites)+length(network.actions));

fixed     = zeros(1,length(network.metabolites)+length(network.actions));
invisible = zeros(1,length(network.metabolites)+length(network.actions));

if exist('sbtab_version','file')

  if ischar(layout_file)
    sbtab_table = sbtab_table_load(layout_file);
  else
    % assume that 'layout_file' is a SBtab object
    sbtab_table = layout_file;
  end

  if flag_KEGG_ids
    A{1} = sbtab_table_get_column(sbtab_table,'Element');
    if sbtab_table_has_column(sbtab_table,'Element_Identifiers_kegg_compound');
      kegg_compound = sbtab_table_get_column(sbtab_table,'Element_Identifiers_kegg_compound');
      ind_exist_kegg_compound = find(cellfun('length',kegg_compound));
      A{1}(ind_exist_kegg_compound) = kegg_compound(ind_exist_kegg_compound);
    end
    if sbtab_table_has_column(sbtab_table,'Element_Identifiers_kegg_reaction');
      kegg_reaction = sbtab_table_get_column(sbtab_table,'Element_Identifiers_kegg_reaction');
      ind_exist_kegg_reaction = find(cellfun('length',kegg_reaction));
      A{1}(ind_exist_kegg_reaction) = kegg_reaction(ind_exist_kegg_reaction);
    end
  else
    A{1} = sbtab_table_get_column(sbtab_table,'Element');
  end

  % positions in columns !Position:x and !Position:y (in older files: !PositionX and !PositionY)
  if sbtab_table_has_column(sbtab_table,'Position_x'),
    A{2} = sbtab_table_get_column(sbtab_table,'Position_x',1);
    A{3} = sbtab_table_get_column(sbtab_table,'Position_y',1);
  else
    A{2} = sbtab_table_get_column(sbtab_table,'PositionX',1);
    A{3} = sbtab_table_get_column(sbtab_table,'PositionY',1);
  end
  A{4} = ones(size(A{2}));

  if flag_fixed
    if sbtab_table_has_column(sbtab_table,'IsFixed');
      A{4} = sbtab_table_get_column(sbtab_table,'IsFixed',1);
    end
  end

  if sbtab_table_has_column(sbtab_table,'IsInvisible');
    A{5} = sbtab_table_get_column(sbtab_table,'IsInvisible',1);
  else
    A{5} = [];
  end

else
  warning('SBtab package is not installed; please install it!');
  fid = fopen(layout_file);
  A   = textscan(fid,'%s%f%f','delimiter','\t');
  fclose(fid);
  A{4} = ones(size(A{2}));
  A{5} = [];

end

% ---

switch flag_KEGG_ids

  case 0
    model_metabolite_names = network.metabolites;
    model_reaction_names   = network.actions;

  case 1

    if isfield(network,'metabolite_KEGGID')
      model_metabolite_names = network.metabolite_KEGGID;
    elseif isfield(network,'Identifiers_kegg_compound')
      model_metabolite_names = network.Identifiers_kegg_compound;
    else
      model_metabolite_names = network.Compound_Identifiers_kegg_compound;
    end

    if length(reaction_KEGG_IDs)
      model_reaction_names = reaction_KEGG_IDs;
    else
      model_reaction_names = network.reaction_KEGGID;
    end

end

if length(unique(model_metabolite_names)) < length(model_metabolite_names)
  if flag_KEGG_ids
    warning('non-unique metabolite names');
  else
    sort(model_metabolite_names)
    error('non-unique metabolite names');
  end
end

if length(unique(model_reaction_names)) < length(model_reaction_names)
  if flag_KEGG_ids
    warning('non-unique reaction names');
  else
    error('non-unique reaction names');
  end

end

l1 = label_names(A{1},model_metabolite_names);
l1(isnan(A{2})) = 0;
ok = find(l1);
x(1,l1(ok)) = A{2}(ok)';
x(2,l1(ok)) = A{3}(ok)';
ll = label_names(model_metabolite_names,A{1});
fixed(find(ll)) = A{4}(ll(find(ll)));
if length(A{5})
  invisible(find(ll)) = A{5}(ll(find(ll)));
end

% ---

l1 = label_names(A{1},model_reaction_names);
l1(isnan(A{2})) = 0;
ok = find(l1);

x(1,length(network.metabolites)+l1(ok)) = A{2}(ok)';
x(2,length(network.metabolites)+l1(ok)) = A{3}(ok)';
ll = label_names(model_reaction_names,A{1});
fixed(length(network.metabolites)+find(ll)) = A{4}(ll(find(ll)));
if length(A{5})
  invisible(length(network.metabolites)+find(ll)) = A{5}(ll(find(ll)));
end

if fill_nans
  x = place_missing_elements(x, network.N);
end

x = x - repmat(offsets',1,size(x,2));

if ~isfield(network,'graphics_par')
  network = netgraph_make_graph(network,[],[]);
end

network.graphics_par.x = x;
network.graphics_par.fixed     = fixed;
network.graphics_par.metinvisible = double(invisible(1:length(model_metabolite_names))==1);
network.graphics_par.actinvisible = double(invisible(length(model_metabolite_names)+1:end)==1);


% ------------------------------------------------------------------

function x = place_missing_elements(x, N)

% Place nodes with unknown (nan) positions: in each round, every missing node that has
% placed neighbours (metabolites and reactions connected by the stoichiometric matrix)
% is put at the mean position of these neighbours, plus a small random shift.
% Nodes without placed neighbours are put around the centre of the layout.

[nm,nr] = size(N);
missing = find(any(isnan(x),1));
if isempty(missing), return; end

placed = all(isfinite(x),1);
adjacency = [sparse(nm,nm), sparse(N~=0); sparse(N'~=0), sparse(nr,nr)];
if any(placed)
  extent = max(max(x(:,placed),[],2) - min(x(:,placed),[],2));
  centre = mean(x(:,placed),2);
else
  extent = 0; centre = [0.5; 0.5];
end
if extent == 0, extent = 1; end

% random numbers from a separate stream (the global random number generator is not changed)
stream = RandStream('mt19937ar','Seed',0);

while ~isempty(missing)
  new_positions = nan(2,length(missing));
  for it = 1:length(missing)
    neighbours = find(adjacency(missing(it),:) & placed);
    if ~isempty(neighbours)
      new_positions(:,it) = mean(x(:,neighbours),2) + 0.04 * extent * randn(stream,2,1);
    end
  end
  done = all(isfinite(new_positions),1);
  if ~any(done), break; end
  x(:,missing(done)) = new_positions(:,done);
  placed(missing(done)) = true;
  missing = missing(~done);
end

if ~isempty(missing)
  x(:,missing) = repmat(centre,1,length(missing)) + 0.1 * extent * randn(stream,2,length(missing));
end
