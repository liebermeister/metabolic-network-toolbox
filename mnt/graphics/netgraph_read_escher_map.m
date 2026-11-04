function [network, report] = netgraph_read_escher_map(network, escher_file, layout_file, options)

% NETGRAPH_READ_ESCHER_MAP Network layout from an Escher map
%
% [network, report] = netgraph_read_escher_map(network, escher_file, layout_file, options)
%
% Reads an Escher map (JSON file, https://escher.github.io), matches its metabolites and
% reactions to the elements of the network, and uses the positions of the Escher nodes as
% node positions: metabolites at their metabolite nodes, reactions at their midpoints.
% Elements without a match in the map are placed near their neighbours in the network.
% The layout is stored in network.graphics_par (the network can be drawn right away) and,
% if layout_file is given, saved as an SBtab position table (see netgraph_print_positions).
% In the table, IsFixed = 1 marks the elements whose positions come from the map.
%
% Model elements and map elements are matched in this order (the first match counts):
%   1. IDs given in options.metabolite_ids and options.reaction_ids
%   2. model IDs (network.metabolites, network.actions) and BiGG IDs in the map; IDs are
%      compared without SBML prefixes 'M_' and 'R_', ignoring upper and lower case and
%      punctuation (e.g. 'M_glc__D_c' matches 'glc__D_c')
%   3. BiGG identifiers of the model (fields 'Identifiers_bigg_metabolite' or 'metabolite_BIGGID',
%      'Identifiers_bigg_reaction' or 'reaction_BIGGID', as read from SBtab columns
%      !Identifiers:bigg.metabolite etc.); metabolite identifiers without compartment (e.g. 'glc__D')
%      match the map metabolite in the compartment given by the end of the model ID (e.g. '_c'),
%      or the only one there is
%   4. names (fields 'metabolite_names', 'reaction_names'), ignoring upper and lower case and
%      punctuation; for metabolites in several compartments, as in 3.
% If a metabolite appears several times in the map, the primary node connected to most
% reactions is used. If a reaction appears several times, the copy closest to its metabolites is used.
%
% Inputs:
%   network      network structure (see 'help metabolic_network'), e.g. from sbtab_to_network
%   escher_file  Escher map (JSON file)
%   layout_file  (optional) name of the SBtab position file to be written ([] for none)
%   options      (optional) struct with fields:
%     metabolite_ids       (default {}) BiGG IDs in the map for the metabolites (cell array, one entry
%                          per metabolite; '' for automatic matching)
%     reaction_ids         (default {}) the same for reactions
%     match_by_names       (default 1) use names for matching (step 4)
%     normalise_positions  (default 1) scale positions to a layout of size 1;
%                          0: Escher coordinates (with y axis pointing upwards)
%     hide_unmatched       (default 0) mark elements that are not in the map as invisible
%                          (graphics_par.metinvisible and actinvisible; column IsInvisible in the file)
%     verbose              (default 0) print a summary of the matching
%
% Outputs:
%   network  network with layout in network.graphics_par
%   report   struct with fields
%              metabolites, reactions  tables with columns ModelID, EscherID, MatchedBy, Copies
%                                      (MatchedBy: 'given', 'ID', 'BiGG identifier', 'name', or '')
%              unused_metabolites, unused_reactions  BiGG IDs of map elements without a model element
%
% Example:
%   network = sbtab_to_network('my_model.tsv', struct('one_sbtab_file', 1));
%   [network, report] = netgraph_read_escher_map(network, 'my_map.json', 'my_model_Position.tsv');
%   netgraph_concentrations(network, [], [], 1);
%   % without metabolites that appear several times in the map (e.g. cofactors):
%   netgraph_draw(network, 'omitmetabolites', report.metabolites.ModelID(report.metabolites.Copies > 1));
%
% See also: netgraph_read_positions, netgraph_print_positions, netgraph_make_graph

eval(default('layout_file','[]','options','struct'));

options_default = struct('metabolite_ids', {{}}, 'reaction_ids', {{}}, 'match_by_names', 1, ...
                         'normalise_positions', 1, 'hide_unmatched', 0, 'verbose', 0);
options = join_struct(options_default, options);

[nm, nr] = size(network.N);
map = escher_map_load(escher_file);

% ------------------------------------------------------------------
% match metabolites and reactions

model_met = struct('ids', {column(network.metabolites)}, ...
                   'given', {model_field(network, {}, nm, options.metabolite_ids)}, ...
                   'bigg', {model_field(network, {'Identifiers_bigg_metabolite','metabolite_BIGGID'}, nm)}, ...
                   'names', {model_field(network, {'metabolite_names'}, nm)});
model_act = struct('ids', {column(network.actions)}, ...
                   'given', {model_field(network, {}, nr, options.reaction_ids)}, ...
                   'bigg', {model_field(network, {'Identifiers_bigg_reaction','reaction_BIGGID'}, nr)}, ...
                   'names', {model_field(network, {'reaction_names'}, nr)});
if ~options.match_by_names
  model_met.names = repmat({''}, nm, 1);
  model_act.names = repmat({''}, nr, 1);
end

[met_match, met_how] = match_elements(model_met, map.metabolites, 1);
[act_match, act_how] = match_elements(model_act, map.reactions,   0);

% ------------------------------------------------------------------
% positions: metabolite nodes; for each reaction, the copy closest to its metabolites

x = nan(2, nm + nr);
for i = find(met_match' > 0)
  x(:,i) = map.metabolites.x(:, met_match(i));
end
for j = find(act_match' > 0)
  copies = map.reactions.copies{act_match(j)};
  x_copies = map.reaction_copies.x(:, copies);
  x_met = x(:, network.N(:,j) ~= 0);
  x_met = x_met(:, all(isfinite(x_met),1));
  best = 1;
  if numel(copies) > 1 && ~isempty(x_met)
    centre = mean(x_met, 2);
    [~, best] = min(sum((x_copies - centre).^2, 1));
  end
  x(:, nm+j) = x_copies(:, best);
end

% Escher coordinates (y axis pointing downwards) -> layout coordinates
if options.normalise_positions
  extent = max(map.box(2) - map.box(1), map.box(4) - map.box(3));
  if extent == 0, extent = 1; end
  x = [(x(1,:) - map.box(1)) / extent; (map.box(4) - x(2,:)) / extent];
else
  x = [x(1,:); -x(2,:)];
end

% ------------------------------------------------------------------
% layout: positions from the map are fixed, other elements are placed near their neighbours

found = find(all(isfinite(x),1));
if isempty(found)
  error('mnt:netgraph:noEscherMatch', ['No element of the network was found in the Escher map %s. ' ...
        'Give the map IDs in options.metabolite_ids and options.reaction_ids'], escher_file);
end
elements = [network.metabolites(:); network.actions(:)];
table_positions = sbtab_table_construct(struct('TableType','Position','TableID','Position'), ...
                    {'Position:x','Position:y'}, {x(1,found), x(2,found)});
table_positions = sbtab_table_add_column(table_positions, 'Element', elements(found), 1);
table_positions = sbtab_table_add_column(table_positions, 'IsFixed', ones(numel(found),1), 1);

if isfield(network, 'graphics_par')
  network = netgraph_read_positions(network, table_positions, [0,0], 1, 0);
else
  network = netgraph_make_graph(network, [], [], [], [], table_positions);
end

if options.hide_unmatched
  network.graphics_par.metinvisible = double(met_match == 0);
  network.graphics_par.actinvisible = double(act_match == 0);
end

if ~isempty(layout_file)
  netgraph_print_positions(network, layout_file, [0,0], 'replace file', 0);
end

% ------------------------------------------------------------------
% report

report.metabolites = match_table(model_met.ids, map.metabolites, met_match, met_how);
report.reactions   = match_table(model_act.ids, map.reactions,   act_match, act_how);
report.unused_metabolites = map.metabolites.ids(setdiff(1:numel(map.metabolites.ids), met_match));
report.unused_reactions   = map.reactions.ids(setdiff(1:numel(map.reactions.ids), act_match));

if options.verbose
  print_summary('Metabolites', met_how, numel(report.unused_metabolites));
  print_summary('Reactions',   act_how, numel(report.unused_reactions));
  if ~isempty(layout_file), fprintf('Layout saved to file %s\n', layout_file); end
end


% =================================================================

function values = model_field(network, field_names, n, given)

  % first field (of length n) found in the network, as a cell array of strings
  % ('' for missing values); or the values given (if not empty)

  values = repmat({''}, n, 1);
  if nargin > 3
    if ~isempty(given)
      if numel(given) ~= n || ~iscell(given)
        error('mnt:netgraph:wrongSize', 'options.metabolite_ids and options.reaction_ids must be cell arrays with one entry per metabolite or reaction');
      end
      values = column(given);
    end
  else
    for k = 1:numel(field_names)
      if isfield(network, field_names{k}) && numel(network.(field_names{k})) == n && iscell(network.(field_names{k}))
        values = column(network.(field_names{k}));
        break
      end
    end
  end
  values(~cellfun(@ischar, values)) = {''};
  % identifiers may contain a prefix, e.g. 'bigg.metabolite:glc__D'
  values = regexprep(values, '^.*[:/]', '');


  % -----------------------------------------------------------------

function [match, how] = match_elements(model, escher, is_metabolite)

  % index of the matching element in the map (0 if none) and the kind of match

  n     = numel(model.ids);
  match = zeros(n,1);
  how   = repmat({''}, n, 1);

  escher_norm  = normalise_id(escher.ids);
  escher_names = normalise_name(escher.names);
  if is_metabolite
    escher_base = regexprep(escher_norm, '_[a-z][a-z0-9]?$', '');   % without compartment
  else
    escher_base = escher_norm;
  end
  model_norm = normalise_id(model.ids);

  % compartments: endings of map IDs (e.g. 'c' in 'glc_d_c'), and of model IDs (if they occur in the map)
  escher_comp = cellfun(@(a,b) regexprep(a(numel(b)+1:end), '^_', ''), escher_norm, escher_base, 'UniformOutput', false);
  model_comp  = repmat({''}, n, 1);
  if is_metabolite
    for c = setdiff(unique(escher_comp(:))', {''})
      ends_with_c = cellfun(@(m) numel(m) > numel(c{1})+1 && strcmp(m(end-numel(c{1}):end), ['_' c{1}]), model_norm);
      model_comp(ends_with_c) = c;
    end
  end

  steps = {'given', normalise_id(model.given), 'exact'; ...
           'ID', model_norm, 'exact'; ...
           'BiGG identifier', normalise_id(model.bigg), 'base'; ...
           'ID', model_norm, 'base'; ...
           'name', normalise_name(model.names), 'name'};

  for s = 1:size(steps,1)
    keys = steps{s,2};
    for i = find(match' == 0)
      key = keys{i};
      if isempty(key), continue; end
      switch steps{s,3}
        case 'exact', candidates = find(strcmp(escher_norm, key));
        case 'base',  candidates = find(strcmp(escher_base, key));
        case 'name',  candidates = find(strcmp(escher_names, key));
      end
      if is_metabolite && ~strcmp(steps{s,3}, 'exact') && ~isempty(model_comp{i})
        % the model ID ends with a compartment: the map element must be in the same compartment
        candidates = candidates(strcmp(escher_comp(candidates), model_comp{i}));
      end
      if numel(candidates) == 1
        match(i) = candidates;
        how{i}   = steps{s,1};
      end
    end
  end


  % -----------------------------------------------------------------

function ids = normalise_id(ids)

  % lower case, without SBML prefixes M_ and R_, compartments '[c]' -> '_c', punctuation -> '_'

  ids = lower(ids);
  ids = regexprep(ids, '^(m|r)_', '');
  ids = regexprep(ids, '\[(\w+)\]$', '_$1');
  ids = regexprep(ids, '[^a-z0-9]+', '_');
  ids = regexprep(ids, '^_|_$', '');


  % -----------------------------------------------------------------

function names = normalise_name(names)

  names = regexprep(lower(names), '[^a-z0-9]', '');


  % -----------------------------------------------------------------

function t = match_table(model_ids, escher, match, how)

  escher_id = repmat({''}, numel(match), 1);
  copies    = zeros(numel(match), 1);
  escher_id(match > 0) = escher.ids(match(match > 0));
  copies(match > 0)    = escher.n_copies(match(match > 0));
  t = table(model_ids, escher_id, how, copies, 'VariableNames', {'ModelID','EscherID','MatchedBy','Copies'});


  % -----------------------------------------------------------------

function print_summary(title, how, n_unused)

  kinds = {'given','ID','BiGG identifier','name'};
  counts = cellfun(@(k) sum(strcmp(how, k)), kinds);
  details = strjoin(arrayfun(@(k) sprintf('%s %d', kinds{k}, counts(k)), find(counts), 'UniformOutput', false), ', ');
  fprintf('%s: %d of %d matched (%s); %d not found in the map; %d map elements not in the model\n', ...
          title, sum(counts), numel(how), details, sum(cellfun(@isempty, how)), n_unused);


  % -----------------------------------------------------------------

function map = escher_map_load(escher_file)

  % metabolites and reactions of an Escher map, with positions
  %   map.metabolites:  ids, names, x (2 x n; position of the chosen node), n_copies
  %   map.reactions:    ids, names, copies (indices in map.reaction_copies), n_copies
  %   map.reaction_copies.x  positions of all reaction copies (midmarkers)
  %   map.box           [xmin xmax ymin ymax] of all nodes

  if ~exist(escher_file, 'file')
    error('mnt:netgraph:fileNotFound', 'Escher map %s not found', escher_file);
  end
  data = jsondecode(fileread(escher_file));
  if iscell(data), data = data{end}; else, data = data(end); end
  if ~isstruct(data) || ~isfield(data, 'nodes') || ~isfield(data, 'reactions')
    error('mnt:netgraph:escherFormat', 'File %s is not an Escher map (no nodes and reactions found)', escher_file);
  end

  % nodes (field names are the node IDs, made valid Matlab names by jsondecode)
  node_keys  = fieldnames(data.nodes);
  nodes      = struct2cell(data.nodes);
  node_type  = cellfun(@(s) s.node_type, nodes, 'UniformOutput', false);
  node_x     = [cellfun(@(s) s.x, nodes), cellfun(@(s) s.y, nodes)]';
  node_bigg  = cellfun(@(s) get_field(s, 'bigg_id', ''), nodes, 'UniformOutput', false);
  node_name  = cellfun(@(s) get_field(s, 'name', ''), nodes, 'UniformOutput', false);
  node_primary = cellfun(@(s) isequal(get_field(s, 'node_is_primary', false), true), nodes);
  node_index = containers.Map(node_keys, num2cell(1:numel(node_keys)));
  node_degree = zeros(numel(nodes),1);   % number of reaction segments at the node

  % reactions: position of the midmarker of each copy
  reactions = struct2cell(data.reactions);
  n_copies  = numel(reactions);
  rxn_bigg  = cell(n_copies,1);
  rxn_name  = cell(n_copies,1);
  rxn_x     = nan(2, n_copies);
  for k = 1:n_copies
    r = reactions{k};
    rxn_bigg{k} = get_field(r, 'bigg_id', '');
    rxn_name{k} = get_field(r, 'name', '');
    ends = [];
    if isfield(r, 'segments') && isstruct(r.segments)
      segments = struct2cell(r.segments);
      for q = 1:numel(segments)
        for node_id = {segments{q}.from_node_id, segments{q}.to_node_id}
          key = matlab.lang.makeValidName(char(node_id{1}));
          if isKey(node_index, key), ends(end+1) = node_index(key); end %#ok<AGROW>
        end
      end
    end
    node_degree(ends) = node_degree(ends) + 1;
    mid = ends(strcmp(node_type(ends), 'midmarker'));
    if ~isempty(mid)
      rxn_x(:,k) = node_x(:, mid(1));
    elseif ~isempty(ends)
      rxn_x(:,k) = mean(node_x(:, ends), 2);
    elseif isfield(r, 'label_x')
      rxn_x(:,k) = [r.label_x; r.label_y];
    end
  end
  has_id = ~cellfun(@isempty, rxn_bigg) & all(isfinite(rxn_x),1)';
  [ids, first, which] = unique(rxn_bigg(has_id), 'stable');
  copy_index = find(has_id);
  map.reactions.ids      = ids;
  map.reactions.names    = rxn_name(copy_index(first));
  map.reactions.copies   = arrayfun(@(u) copy_index(which == u)', (1:numel(ids))', 'UniformOutput', false);
  map.reactions.n_copies = cellfun(@numel, map.reactions.copies);
  map.reaction_copies.x  = rxn_x;

  % metabolites: of several nodes, the primary node with most segments
  is_met = find(strcmp(node_type, 'metabolite') & ~cellfun(@isempty, node_bigg));
  [ids, ~, which] = unique(node_bigg(is_met), 'stable');
  map.metabolites.ids      = ids;
  map.metabolites.names    = cell(numel(ids),1);
  map.metabolites.x        = nan(2, numel(ids));
  map.metabolites.n_copies = zeros(numel(ids),1);
  for u = 1:numel(ids)
    copies = is_met(which == u);
    score  = 1000 * node_primary(copies) + node_degree(copies);
    [~, best] = max(score);
    map.metabolites.names{u}    = node_name{copies(best)};
    map.metabolites.x(:,u)      = node_x(:, copies(best));
    map.metabolites.n_copies(u) = numel(copies);
  end

  map.box = [min(node_x(1,:)), max(node_x(1,:)), min(node_x(2,:)), max(node_x(2,:))];


  % -----------------------------------------------------------------

function value = get_field(s, name, default_value)

  if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
  else
    value = default_value;
  end
