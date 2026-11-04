function [names, positions] = netgraph_print_positions(network, layout_file, offsets,policy,flag_KEGG_ids,flag_element_names,flag_fixed,flag_invisible)

% NETGRAPH_PRINT_POSITIONS Save a network layout to a file
%
% [names, positions] = netgraph_print_positions(network,layout_file,offsets,policy,flag_KEGG_ids,flag_element_names,flag_fixed,flag_invisible)
%
% Argument 'policy' (used if the file exists already):
%   'replace elements' (default)  write all model elements; keep file elements that are not in the model
%   'add nonexistent'             keep all file elements; add model elements that are not in the file
%   'replace file'                write the model elements only
%
% See also: netgraph_read_positions

eval(default('layout_file','[]','offsets','[]','policy','[]','flag_KEGG_ids','1','flag_element_names','1','flag_fixed','1','flag_invisible','1'));

if flag_KEGG_ids
  if ~[isfield(network,'metabolite_KEGGID') * isfield(network,'reaction_KEGGID')]
    %warning('KEGG IDs missing in model - using model element IDs in position file');
    flag_KEGG_ids=0;
  end
end


offsets = column(offsets)';
if isempty(offsets), offsets = [0 0]; end
if isempty(policy), policy = 'replace elements'; end

[nr,nm] = network_numbers(network);

metnames = network.metabolites;
actnames = network.actions;

if flag_KEGG_ids
  if ~[isfield(network,'metabolite_KEGGID') * isfield(network,'reaction_KEGGID')]
    error('KEGG IDs missing in model');
  end
  kegg_metnames = network.metabolite_KEGGID;
  kegg_actnames = network.reaction_KEGGID;
else
  kegg_metnames = repmat({''},nm,1);
  kegg_actnames = repmat({''},nr,1);
end

metpositions = network.graphics_par.x(:,1:nm) + repmat(offsets',1,nm);
actpositions = network.graphics_par.x(:,nm+1:end) + repmat(offsets',1,nr);

names         = [metnames; actnames];
kegg_compound = [column(kegg_metnames); repmat({''},nr,1)];
kegg_reaction = [repmat({''},nm,1); column(kegg_actnames)];
positions     = [metpositions, actpositions];
if isfield(network.graphics_par,'fixed')
  fixed       = column(network.graphics_par.fixed);
else
  fixed       = zeros(nm+nr,1);
end
if isfield(network.graphics_par,'metinvisible')
  invisible   = [column(network.graphics_par.metinvisible); column(network.graphics_par.actinvisible)];
else
  invisible   = zeros(nm+nr,1);
end

% -------------------------------------------------

if length(layout_file)

  if ~exist(layout_file,'file'), policy = 'replace file'; end

  switch policy

    case 'replace file'

    case {'add nonexistent', 'replace elements'}

      t             = sbtab_table_load(layout_file);
      elements      = sbtab_table_get_column(t,'Element');
      % positions in columns !Position:x and !Position:y (in older files: !PositionX and !PositionY)
      if sbtab_table_has_column(t,'Position_x'),
        x           = cell_string2num(sbtab_table_get_column(t,'Position_x'));
        y           = cell_string2num(sbtab_table_get_column(t,'Position_y'));
      else
        x           = cell_string2num(sbtab_table_get_column(t,'PositionX'));
        y           = cell_string2num(sbtab_table_get_column(t,'PositionY'));
      end

      kegg_compound_in = repmat({''},length(elements),1);
      kegg_reaction_in = repmat({''},length(elements),1);

      if flag_KEGG_ids
        if sbtab_table_has_column(t,'Element_Identifiers_kegg_compound')
          kegg_compound_in = sbtab_table_get_column(t,'Element_Identifiers_kegg_compound');
        end
        if  sbtab_table_has_column(t,'Element_Identifiers_kegg_reaction')
          kegg_reaction_in = sbtab_table_get_column(t,'Element_Identifiers_kegg_reaction');
        end
      end

      fixed_in  = zeros(length(elements),1);
      if flag_fixed
        if sbtab_table_has_column(t,'IsFixed')
          fixed_in = column(sbtab_table_get_column(t,'IsFixed',1));
        end
      end

      invisible_in  = zeros(length(elements),1);
      if flag_invisible
        if  sbtab_table_has_column(t,'IsInvisible')
          invisible_in = column(sbtab_table_get_column(t,'IsInvisible',1));
        end
      end

      x = column(x); y = column(y);

      switch policy

        case 'add nonexistent'
          % keep all elements from the file; add model elements that are missing in the file
          ll            = label_names(names,elements);
          ind_keep      = find(ll==0);
          names         = [elements;         names(ind_keep)];
          kegg_compound = [kegg_compound_in; kegg_compound(ind_keep)];
          kegg_reaction = [kegg_reaction_in; kegg_reaction(ind_keep)];
          positions     = [[x';y'],          positions(:,ind_keep) ];
          fixed         = [fixed_in;         fixed(ind_keep)];
          invisible     = [invisible_in;     invisible(ind_keep)];

        case 'replace elements'
          % keep file elements that are missing in the model; write all model elements
          ll            = label_names(elements,names);
          ind_keep      = find(ll==0);
          names         = [elements(ind_keep);         names];
          kegg_compound = [kegg_compound_in(ind_keep); kegg_compound];
          kegg_reaction = [kegg_reaction_in(ind_keep); kegg_reaction];
          positions     = [[x(ind_keep)';y(ind_keep)'], positions];
          fixed         = [fixed_in(ind_keep);         fixed];
          invisible     = [invisible_in(ind_keep);     invisible];
      end

  end

  % ----------------------------------------------------------------

  sbtab_out = sbtab_table_construct( struct('TableType','Position','TableID','Position'), ...
                      {'Position:x','Position:y'},  {positions(1,:),positions(2,:)});

  if flag_element_names
    sbtab_out = sbtab_table_add_column(sbtab_out,'Element',names,1);
  end

  if flag_KEGG_ids
    sbtab_out = sbtab_table_add_column(sbtab_out,'Element:Identifiers:kegg:compound',kegg_compound,1);
    sbtab_out = sbtab_table_add_column(sbtab_out,'Element:Identifiers:kegg:reaction',kegg_reaction,1);
  end

  if flag_fixed
    sbtab_out = sbtab_table_add_column(sbtab_out,'IsFixed',fixed,1);
  end

  if flag_invisible
    sbtab_out = sbtab_table_add_column(sbtab_out,'IsInvisible',invisible,1);
  end

  if length(layout_file)
    sbtab_table_save(sbtab_out, struct('filename', layout_file));
  end

end
