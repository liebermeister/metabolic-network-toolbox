% NETGRAPH_DRAW Draw a network
%
% network = netgraph_draw(network, gp)
% network = netgraph_draw(network, 'name1', value1, 'name2', value2, ...)
%
% Draw a metabolic network according to settings in the field 'graphics_par'.
% These settings can be overridden by graphics options, given as a struct options
% (e.g. struct('metvalues', c, 'colorbar', 1)), as name-value pairs, or as a struct
% followed by name-value pairs. All graphics functions in this folder take options
% in the same way and pass them on to netgraph_draw; the options are listed below.
%
% Option names are checked: names that differ only in upper and lower case are corrected
% (e.g. 'fontsize' -> 'FontSize'); unknown names give a warning (once per Matlab session),
% which can be switched off by warning('off','mnt:netgraph:unknownOption')
%
% Example: netgraph_draw(network, 'metvalues', c, 'arrowstyle', 'fluxes', 'arrowvalues', v)
%
%PARAMETER SETTINGS:
%
%METABOLITE SYMBOLS
% metvalues
% metvalues_std
% metstyle       {'fixed','box','box_std','none'}
% metcolstyle    {'fixed','values'}
% metvaluesmax
% metvaluesmin
% metprintnames
% metprintvalues
% metcolors       (optional, overrides all other color specifications)
% metinvisible    (bit vector)
% omitmetabolites (list of metabolite names) metabolites to be omitted in graphics
% squaresize ( also holds for reaction symbols)
%                (all sizes and offsets refer to a layout of size 1; see 'scale_sizes')
% no_points ( also holds for reaction symbols)
% flag_edges (holds also for reaction and regulation symbols)
% circle_shift [0,0] additional shift of circle symbols
%
%REACTION SYMBOLS
% actvalues
% actvalues_std
% actstyle       {'fixed','box','box_std','none'}
% actcolstyle    {'fixed','values'}
% actcolors      (optional, overrides all other color specifications)
% actinvisible    (bit vector)
% actvaluesmax
% actvaluesmin
% actprintnames
% actprintvalues
% omitreactions   (list of reaction names) reactions to be omitted in graphics
%
%STYLE OF ARROWS
% arrowstyle       {'fluxes','directions','none'}
% arrowvalues      a value of 1 will be shown in standard size (given by arrowsize);
%                  unknown values (nan) are shown by grey diamonds
% arrowvaluesmax   threshold value; any larger values will be replaced by this value
% arrowvaluesmin   threshold value; any smaller (or more negative) values will be replaced by this value
% arrow_stoichiometries [0,1] scale arrows with stoichiometric coefficients? (default 0)
% arrowsize
% arrow_shift      [0..1]
% arrowprintvalues (overrides actprintvalues);
% straightlines    (default 1) 0: draw curved lines that pass smoothly through each reaction
% arrowcolor
% flag_triangle_edges (default 0)
% single_arrow     [0,1] Draw only one arrow per reaction, on top of the reaction symbol
%
%STYLE OF EDGE SYMBOLS
%
% edgevalues     (sparse matrix, size of stoichiometric matrix N)
% edgevalues_std (sparse matrix, size of stoichiometric matrix N)
% edgestyle      {'normal',fixed'}; 'fixed' plot all symbols in the same size
% edgevaluesmax
% edgevaluesmin
% suppress lines flag - if set, no lines are drawn
%
%REGULATION ARROWS
% show_regulation (flag) - requires options entry "regulation_matrix" to be given
% show_regulationvalues (flag)
% regulationvalues (sparse matrix of size of regulation_matrix)
% regulationstyle {'normal',fixed'} 'fixed' plot all symbols in the same size
% rlinecolor
%
%TEXT ITEMS
% chop_names     (length)
% split_names (length)
% FontSize
% fontangle       rotation angle of all text items (can also be given as 'Rotation')
% text_offset     [dx,dy] position of name labels relative to the nodes
%                 (default: to the right of the symbols, vertically centred)
% value_text_offset
%
%COLORSCALE AND LINES
% colorbar (0/1)
% colorbar_numbers
% colorbar_fontsize
% colorbar_location (string, e.g. 'South' (default), 'East', 'West'); colorbars are placed outside the network
% linecolor
% linewidth
% colormap        (default rb_colors)
% black_and_white (default 0)
% showsign        (default 1) 1: colours show values on a symmetric scale [-M,M], with
%                 M = max(|valuesmin|,|valuesmax|); symbol sizes show absolute values
%                 0: colours and sizes show values on the scale [valuesmin, valuesmax]
% shade_long_lines  (value) if ~=0 : length above which  lines are shaded
% shade_cofactor_lines  (bitstring denoting the cofactors)
%
% canvas_position [xmin ymin xsize ysize] values for figure "Position" property
% canvas (flag, default []): show light grey background; or rgb vector with background color
% figure_axis: axis limits [xmin xmax ymin ymax] (default: limits that contain all nodes and labels)
% figure_position
% subplot_position
% YDir ('normal','reverse')
% hold_on (flag)
% keep_subplot (flag, default 0): do not change the position of the axes.
%     Without this flag, a new axes object (e.g. after 'figure' or 'clf') is enlarged to fill the figure;
%     axes created by 'subplot' or placed by the user keep their position in any case
% scale_sizes (flag, default 1): sizes (squaresize, arrowsize, text offsets, circle_shift, shade_long_lines)
%     are multiplied by the size of the layout (largest extent of the node positions), so the
%     graphics look similar for layouts in any units. Set to 0 to use sizes in layout units
%
% See also: netgraph_make_graph, netgraph_concentrations, netgraph_move

function network = netgraph_draw(network,varargin)

if sum(size(network.N)) > 10000
  warning('mnt:netgraph:tooLarge', 'Network too large for drawing'); return
end

p = netgraph_draw_set_parameters(network,varargin);

n_met = length(network.metabolites);
n_act = length(network.actions);
x     = p.x;

% sizes are given for a layout of size 1; scale them to the size of the actual layout
size_fields = {'squaresize','arrowsize','text_offset','value_text_offset','circle_shift','shade_long_lines'};
unscaled = struct;
for it = 1:length(size_fields)
  unscaled.(size_fields{it}) = p.(size_fields{it});
  p.(size_fields{it})        = p.layout_size * p.(size_fields{it});
end

if ~isempty(p.figure_position), set(gcf,'Position',p.figure_position); end
if ~isempty(p.subplot_position), subplot('position',p.subplot_position); end

ax = gca;
had_content = p.hold_on==1 && ~isempty(get(ax,'Children'));
if p.hold_on==1
  hold(ax,'on');
else
  cla(ax); colorbar(ax,'off');
  set(ax, 'XScale', 'linear', 'YScale', 'linear');   % (e.g. after a plot with log scale)
end

% a new axes object (at the default position) is enlarged to fill the figure;
% axes created by 'subplot' or placed by the user keep their position
if ~p.keep_subplot && ~had_content && isequal(round(get(ax,'Position'),4), round(get(gcf,'DefaultAxesPosition'),4))
  set(ax,'Position',[0.02 0.03 0.96 0.9]);
end

hold(ax,'on');

if ~isempty(p.canvas)
  if length(p.canvas) ~= 3, p.canvas = [.98 .98 .98]; end
  set(gcf,'Color',p.canvas);
else
  set(gcf,'Color',[1 1 1]);
end

% Graphics elements of one kind are drawn as one graphics object (with tags such as
% 'netgraph_edges' or 'netgraph_metabolites'), because single objects for each element
% make drawing slow in large networks

draw_arrows = p.show_arrowvalues && any(strcmp(p.arrowstyle, {'fluxes','directions'}));
nan_size    = 0.4 * p.arrowsize;   % size of the symbols for unknown fluxes

% ------------------------------------------------------------------------------
% edges between metabolites and reactions: points along each edge, from the metabolite to the reaction

[edge_met, edge_act] = find(p.m(p.metindshow, n_met + p.actindshow));
edge_met = reshape(p.metindshow(edge_met), [], 1);
edge_act = reshape(p.actindshow(edge_act), [], 1);
n_edges  = length(edge_met);
stoich   = full(p.N(sub2ind(size(p.N), edge_met, edge_act)));

edge_lines = cell(n_edges,1);
edge_mid   = zeros(2,n_edges);   % midpoints of the edges (positions of edge symbols)
for k = 1:n_edges
  x_met = x(:,edge_met(k));
  x_act = x(:,n_met+edge_act(k));
  if p.straightlines
    edge_lines{k} = [x_met, x_act];
  else
    % curved lines: each edge is a circular arc that continues towards the other side of the reaction
    % (substrate edges bend towards the products, product edges towards the substrates)
    if stoich(k) < 0
      other_side = find(p.N(:,edge_act(k)) > 0);
    else
      other_side = find(p.N(:,edge_act(k)) < 0);
    end
    height = 0;
    if ~isempty(other_side)
      height = arc_height(x_met, x_act, mean(x(:,other_side),2));
    end
    edge_lines{k} = fliplr(arc(x_met, x_act, height, 0.5));
  end
  edge_mid(:,k) = point_on_line(edge_lines{k}, 0.5);
end

% ------------------------------------------------------------------------------
% lines (long lines and lines to cofactors can be shaded)

if ~p.suppress_lines
  shaded = false(n_edges,1);
  if p.shade_long_lines
    shaded = sum((x(:,edge_met) - x(:,n_met+edge_act)).^2, 1)' > p.shade_long_lines^2;
  end
  if ~isempty(p.shade_cofactor_lines)
    shaded = shaded | reshape(p.shade_cofactor_lines(p.metabolite_mapping(edge_met)) ~= 0, [], 1);
  end
  draw_lines(ax, edge_lines(~shaded), p.linecolor, p.linewidth, 'netgraph_edges');
  draw_lines(ax, edge_lines(shaded), 0.2 * p.linecolor + 0.8 * [1 1 1], 1, 'netgraph_shaded_edges');
end

% little dots at all nodes
if ~p.no_points
  dd = [p.metindshow, n_met + p.actindshow];
  plot(ax, x(1,dd), x(2,dd), '.', 'Color', p.linecolor, 'MarkerSize', 6, 'Tag', 'netgraph_nodes');
end

% ------------------------------------------------------------------------------
% arrows on the edges (forward = 1: from metabolite to reaction; -1: opposite direction)

if draw_arrows && ~p.single_arrow && n_edges
  flux    = reshape(p.norm_arrowvalues(edge_act), [], 1);
  forward = -sign(stoich) .* sign(reshape(p.arrowvalues(edge_act), [], 1));
  len     = p.arrowsize * sqrt(abs(flux));
  if p.arrow_stoichiometries, len = len .* abs(stoich); end
  shift   = p.arrow_shift * (forward == 1) + (1 - p.arrow_shift) * (forward == -1) + 0.5 * isnan(forward);
  pos = zeros(2,n_edges);
  phi = zeros(n_edges,1);
  for k = 1:n_edges
    [pos(:,k), direction] = point_on_line(edge_lines{k}, shift(k));
    phi(k) = atan2(direction(2), direction(1));
  end
  draw_arrow_symbols(ax, pos, phi, len, forward, p.arrowcolor, p.flag_triangle_edges, nan_size);
end

% ---------------------------------------------------
% edge values

if p.show_edgevalues && n_edges
  has_value = full(network.N(sub2ind(size(network.N), edge_met, edge_act))) ~= 0;
  ind       = sub2ind(size(p.norm_edgevalues), edge_met(has_value), edge_act(has_value));
  values    = full(p.norm_edgevalues(ind));
  colors    = get_color(values, -1, 1, p.colormap);
  if strcmp(p.edgestyle, 'fixed'), values(isfinite(values)) = 1; end
  draw_symbols(ax, 'diamond', edge_mid(:,has_value), abs(values), full(p.norm_edgevalues_std(ind)), ...
               p.squaresize, colors, p.flag_edges, 'netgraph_edge_symbols');
  if p.edgeprintvalues
    offset = 0.5 * p.squaresize;
    text(ax, edge_mid(1,has_value) + offset, edge_mid(2,has_value) + offset, value_strings(full(p.edgevalues(ind))), ...
         'FontSize', p.FontSize, 'Rotation', p.fontangle);
  end
end

% ---------------------------------------------------
% regulation arrows (arcs from the metabolite to the reaction)

if p.show_regulation
  [ind_r, ind_m] = find(network.regulation_matrix);
  keep  = ismember(ind_r, p.actindshow) & ismember(ind_m, p.metindshow);
  ind_r = ind_r(keep);
  ind_m = ind_m(keep);
  n_reg = length(ind_r);
  reg_lines = cell(n_reg,1);
  reg_pos   = zeros(2,n_reg);
  for k = 1:n_reg
    [reg_lines{k}, reg_pos(1,k), reg_pos(2,k)] = arc(x(:,n_met+ind_r(k)), x(:,ind_m(k)), 0.2, 0.5);
  end
  ind = sub2ind(size(network.regulation_matrix), ind_r, ind_m);
  activation = full(network.regulation_matrix(ind)) > 0;
  if ~p.suppress_lines
    draw_lines(ax, reg_lines(activation),  [0 0.7 1], 0.5, 'netgraph_activation');
    draw_lines(ax, reg_lines(~activation), [1 0 0],   0.5, 'netgraph_inhibition');
  end
  if p.show_regulationvalues
    ind    = sub2ind(size(p.norm_regulationvalues), ind_r, ind_m);
    values = full(p.norm_regulationvalues(ind));
    colors = get_color(values, -1, 1, p.colormap);
    if strcmp(p.regulationstyle, 'fixed'), values(isfinite(values)) = 1; end
    draw_symbols(ax, 'diamond', reg_pos, abs(values), full(p.norm_regulationvalues_std(ind)), ...
                 p.squaresize, colors, p.flag_edges, 'netgraph_regulation_symbols');
  end
  if p.regulationprintvalues
    offset = 0.5 * p.squaresize;
    text(ax, reg_pos(1,:) + offset, reg_pos(2,:) + offset, value_strings(full(p.regulationvalues(ind))), ...
         'FontSize', p.FontSize, 'Rotation', p.fontangle);
  end
end

% ---------------------------------------------------
% metabolites

if p.show_metvalues
  switch p.metcolstyle
    case 'fixed',  c = repmat(p.metcol,n_met,1);
    case 'values', c = get_color(p.norm_metvalues,-1,1,p.colormap);
  end
  if ~isempty(p.metcolors)
    c = p.metcolors(p.metabolite_mapping,:);
  end
  values = p.norm_metvalues(:);
  stds   = p.norm_metvalues_std(:);
  if p.showsign == 0
    values = (values + 1) / 2;
    stds   = 0.5 * stds;
  end
  values = abs(values);
  switch p.metstyle
    case 'box',   stds = 0 * stds;
    case 'fixed', values(isfinite(values)) = 1;
  end
  ind   = p.metindshow;
  known = isfinite(values(ind))';
  if p.suppress_lines, ind = ind(known); known = known(known); end
  pos = x(:,ind);
  pos(:,known) = pos(:,known) + p.circle_shift(:);   % symbols for unknown values are not shifted
  draw_symbols(ax, 'circle', pos, values(ind), stds(ind), p.squaresize, c(ind,:), p.flag_edges, 'netgraph_metabolites');
end

% ---------------------------------------------------
% reactions

if p.show_actvalues
  switch p.actcolstyle
    case 'fixed',  c = repmat(p.actcol,n_act,1);
    case 'values', c = get_color(p.norm_actvalues,-1,1,p.colormap);
  end
  if ~isempty(p.actcolors)
    c = p.actcolors(p.reaction_mapping,:);
  end
  values = abs(p.norm_actvalues(:));
  stds   = p.norm_actvalues_std(:);
  switch p.actstyle
    case 'box',   stds = 0 * stds;
    case 'fixed', values(isfinite(values)) = 1;
  end
  ind   = p.actindshow;
  known = isfinite(values(ind))';
  pos   = x(:,n_met+ind);
  pos(:,known) = pos(:,known) + p.circle_shift(:);
  draw_symbols(ax, 'square', pos, values(ind), stds(ind), p.squaresize, c(ind,:), p.flag_edges, 'netgraph_reactions');
end

% ---------------------------------------------------
% single arrows: one arrow per reaction, on the edge from an incoming metabolite
% (a substrate for positive fluxes, a product for negative fluxes);
% for unknown fluxes (nan), a grey diamond is shown on the reaction

if draw_arrows && p.single_arrow
  pos = zeros(2,0);
  phi = zeros(1,0);
  len = zeros(1,0);
  for ind_act = p.actindshow
    flux = p.norm_arrowvalues(ind_act);
    if flux == 0, continue; end
    if flux > 0 || isnan(flux)
      incoming = find(p.N(:,ind_act) < 0);
    else
      incoming = find(p.N(:,ind_act) > 0);
    end
    incoming = intersect(incoming, p.metindshow);
    if isempty(incoming), continue; end
    x_met = x(:,incoming(1));
    x_act = x(:,n_met+ind_act);
    l = p.arrowsize * sqrt(abs(flux));
    if p.arrow_stoichiometries, l = l * abs(p.N(incoming(1),ind_act)); end
    pos(:,end+1) = x_met + p.arrow_shift * (x_act - x_met);
    phi(end+1)   = atan2(x_act(2) - x_met(2), x_act(1) - x_met(1));
    len(end+1)   = l;
  end
  draw_arrow_symbols(ax, pos, phi, len, ones(size(len)), p.arrowcolor, p.flag_triangle_edges, nan_size);
end

% ---------------------------------------------------
% text

text_style = {'FontSize', p.FontSize, 'Rotation', p.fontangle};

if p.metprintnames
  ind = p.metindshow;
  text(ax, x(1,ind) + p.text_offset(1), x(2,ind) + p.text_offset(2), p.metnames(ind), text_style{:});
end

if p.actprintnames
  ind = p.actindshow;
  text(ax, x(1,n_met+ind) + p.text_offset(1), x(2,n_met+ind) + p.text_offset(2), p.actnames(ind), text_style{:});
end

offset = 0.5 * p.squaresize * [1,1] + p.value_text_offset;

if p.metprintvalues
  if isempty(p.metvalues), error('mnt:netgraph:missingValues', 'No metabolite values given'); end
  strings = value_strings(p.metvalues, p.metvalues_std);
  ind = p.metindshow;
  text(ax, x(1,ind) + offset(1), x(2,ind) + offset(2), strings(ind), text_style{:});
end

if p.actprintvalues
  if isempty(p.actvalues), error('mnt:netgraph:missingValues', 'No reaction values given'); end
  strings = value_strings(p.actvalues, p.actvalues_std);
  ind = p.actindshow;
  text(ax, x(1,n_met+ind) + offset(1), x(2,n_met+ind) + offset(2), strings(ind), text_style{:});
end

if p.arrowprintvalues
  if isempty(p.arrowvalues), error('mnt:netgraph:missingValues', 'No arrow values given'); end
  strings = value_strings(p.arrowvalues, p.arrowvalues_std);
  ind = p.actindshow;
  text(ax, x(1,n_met+ind) + offset(1), x(2,n_met+ind) + offset(2), strings(ind), text_style{:});
end

% ---------------------------------------------------
% general layout

set(ax,'YDir',p.YDir);
axis(ax,'off');
daspect(ax,[1 1 1]);

if ~isempty(p.figure_axis)
  axis(ax,p.figure_axis);
else
  set_axis_limits(ax, x(:,[p.metindshow, n_met+p.actindshow]), p.squaresize, had_content);
end

if p.colorbar
  if isempty(p.colorbar_numbers)
    if ~isempty(p.metvaluesmax)
      p.colorbar_numbers = [p.metvaluesmin,p.metvaluesmax];
    elseif ~isempty(p.actvaluesmax)
      p.colorbar_numbers = [p.actvaluesmin,p.actvaluesmax];
    elseif p.edgevaluesmax
      p.colorbar_numbers = [p.edgevaluesmin,p.edgevaluesmax];
    end
  end
  if max(p.colorbar_numbers) > min(p.colorbar_numbers)
    this_colorbar(ax, p.colorbar_numbers, p.colormap, p.colorbar_location, p.colorbar_fontsize);
  end
end

% return the size parameters as given (not scaled)
for it = 1:length(size_fields)
  p.(size_fields{it}) = unscaled.(size_fields{it});
end

network.graphics_par = p;

if ~isempty(p.canvas_position)
  set(gcf,'Position',p.canvas_position);
end


% =================================================================
% Drawing functions

function draw_lines(ax, lines, color, width, tag)

  % polygonal lines (cell array of 2 x n matrices), drawn as one line object

  if isempty(lines), return; end
  xy = cellfun(@(l) [l, nan(2,1)], lines(:)', 'UniformOutput', false);
  xy = [xy{:}];
  line(ax, xy(1,:), xy(2,:), 'Color', color, 'LineWidth', width, 'Tag', tag);


  % -----------------------------------------------------------------

function draw_polygons(ax, outline, pos, sizes, colors, edge_color, tag)

  % copies of a polygon (outline: 2 x k vertices, for size 1) at positions pos (2 x n),
  % scaled by sizes (n values), with face colors (n x 3), drawn as one patch object

  if isempty(sizes), return; end
  X = pos(1,:) + outline(1,:)' * sizes(:)';
  Y = pos(2,:) + outline(2,:)' * sizes(:)';
  patch(ax, 'XData', X, 'YData', Y, 'FaceVertexCData', colors, 'FaceColor', 'flat', ...
        'EdgeColor', edge_color, 'Tag', tag);


  % -----------------------------------------------------------------

function draw_symbols(ax, shape, pos, values, stds, squaresize, colors, flag_edges, tag)

  % symbols ('circle', 'square', or 'diamond') at positions pos (2 x n), with sizes given by
  % values (between 0 and 1) and colours (n x 3); standard deviations (stds) are shown by a
  % lighter ring around the symbols. Unknown values (nan) are shown by small grey symbols

  switch shape
    case 'circle'
      t = linspace(0, 2*pi, 41); t(end) = [];
      outline = 0.5 * [cos(t); sin(t)];        nan_value = 0.2;
    case 'square'
      outline = 0.5 * [-1 1 1 -1; -1 -1 1 1];  nan_value = 0.2;
    case 'diamond'
      outline = 2/pi * [-1 0 1 0; 0 -1 0 1];   nan_value = 0.5;
  end

  values = values(:);
  stds   = stds(:);
  if isempty(stds), stds = zeros(size(values)); end
  std_colors = 0.7 + 0.3 * colors;

  unknown = isnan(values);
  values(unknown)       = nan_value;
  colors(unknown,:)     = repmat([.7 .7 .7], sum(unknown), 1);
  std_colors(unknown,:) = 0;
  stds(isnan(stds))     = 0.1;

  if flag_edges, edge_color = [.4 .4 .4]; else, edge_color = 'none'; end

  ring = stds ~= 0;
  draw_polygons(ax, outline, pos(:,ring), squaresize * (values(ring) + stds(ring)), std_colors(ring,:), edge_color, [tag '_std']);
  sizes = squaresize * max(0, values);
  shown = sizes > 0;
  draw_polygons(ax, outline, pos(:,shown), sizes(shown), colors(shown,:), edge_color, tag);


  % -----------------------------------------------------------------

function draw_arrow_symbols(ax, pos, phi, len, forward, color, flag_edges, nan_size)

  % arrows (triangles) at positions pos (2 x n) with lengths len, pointing in direction phi
  % (angle; forward = 1) or in the opposite direction (forward = -1);
  % unknown values (nan) are shown by grey diamonds of width 2 * nan_size

  len     = len(:)';
  forward = forward(:)';
  phi     = phi(:)';
  if flag_edges, edge_color = [0 0 0]; else, edge_color = 'none'; end

  unknown = isnan(len) | isnan(forward);
  known   = ~unknown & len > 0 & forward ~= 0;

  if any(known)
    triangle = [-0.5 0.5 -0.5; -0.3 0 0.3];   % pointing in x direction, length 1
    s = forward(known) .* len(known);
    c = cos(phi(known));
    d = sin(phi(known));
    X = pos(1,known) + triangle(1,:)' * (s .* c) - triangle(2,:)' * (s .* d);
    Y = pos(2,known) + triangle(1,:)' * (s .* d) + triangle(2,:)' * (s .* c);
    patch(ax, 'XData', X, 'YData', Y, 'FaceColor', color, 'EdgeColor', edge_color, 'Tag', 'netgraph_arrows');
  end

  if any(unknown)
    diamond = nan_size * [-1 0 1 0; 0 -1 0 1];
    X = pos(1,unknown) + diamond(1,:)';
    Y = pos(2,unknown) + diamond(2,:)';
    patch(ax, 'XData', X, 'YData', Y, 'FaceColor', [.7 .7 .7], 'EdgeColor', edge_color, 'Tag', 'netgraph_unknown_fluxes');
  end


  % -----------------------------------------------------------------

function strings = value_strings(values, values_std)

  % numbers as text labels (empty labels for nan values), with standard deviations if given

  values = full(values(:));
  strings = arrayfun(@(v) num2str(v,3), values, 'UniformOutput', false);
  if nargin > 1 && any(values_std(:) ~= 0)
    stds    = arrayfun(@(v) num2str(v,2), full(values_std(:)), 'UniformOutput', false);
    strings = strcat(strings, {' \pm '}, stds);
  end
  strings(isnan(values)) = {''};


  % -------------------------
  % compute positions on an arc: x3 all points for drawing;
  % arcx,arcy: positions of one point (used for arrow)

function [x3,arcx,arcy] = arc(x1,x2,height,arrow_shift)

  s = 0:0.1:pi;
  x3 = repmat(mean([x1,x2],2),1,length(s)) ...
       + 0.5 * [x2-x1] * cos(s) + height * 0.5 * [0 -1; 1 0] * [x2-x1] * sin(s);
  nn = max(1,ceil(arrow_shift * size(x3,2)));
  arcx = x3(1,nn);
  arcy = x3(2,nn);


  % -----------------------------------------------------------------

function height = arc_height(x1,x2,x3)

  % height parameter for 'arc' (bulge relative to half the chord length) such that
  % the arc from x1 to x2 follows the circle through x1, x2, and x3 (in this order).
  % This arc lies on the side of the chord opposite to x3.
  % Returns 0 if the three points lie on a straight line; values are limited to [-1,1] (half circle)

  R90 = [0 -1; 1 0];
  d   = x2 - x1;
  c   = x3 - 0.5*(x1+x2);
  chord = sqrt(sum(d.^2));
  side  = (R90*d)' * c;   % side of the chord on which x3 lies
  if chord == 0 || abs(side) <= 10^-9 * chord * sqrt(sum(c.^2))
    height = 0; return
  end
  % centre of the circle: (x1+x2)/2 + t * R90 * d
  t = (sum(c.^2) - 0.25*chord^2) / (2*side);
  radius = sqrt(0.25 + t^2) * chord;
  if t * side > 0
    % centre on the same side as x3: the arc is the shorter one
    sagitta = radius - abs(t) * chord;
  else
    sagitta = radius + abs(t) * chord;
  end
  height = - sign(side) * 2 * sagitta / chord;
  height = max(-1, min(1, height));


  % -----------------------------------------------------------------

function [x_point, direction] = point_on_line(x_line, fraction)

  % point at a given fraction of the arc length along a polygonal line,
  % and the direction of the line at this point

  segments = diff(x_line,1,2);
  lengths  = sqrt(sum(segments.^2,1));
  total    = sum(lengths);
  if total == 0
    x_point = x_line(:,1); direction = [1;0]; return
  end
  cum = [0, cumsum(lengths)];
  target = min(max(fraction,0),1) * total;
  k = find(cum(2:end) >= target, 1);
  if isempty(k), k = length(lengths); end
  u = (target - cum(k)) / max(lengths(k),eps);
  x_point   = x_line(:,k) + u * segments(:,k);
  direction = segments(:,k);


  % -------------------------------------------------

function name2 = split_name(name,l)

  name2 = name;
  if length(name) > l
    ll = ceil(length(name)/2);
    name  = [name repmat(' ',1,ll - mod(length(name),ll))];
    name1 = [reshape(name,ll,length(name)/ll)', repmat('\n',length(name)/ll,1)];
    name1 = reshape(name1',prod(size(name1)),1)';
    name2 = sprintf(name1(1:end-2));
  end


  % -------------------------------------------------

function [normvalues, normvalues_std, valuesmax, valuesmin] = normalise(values,values_std,valuesmax,valuesmin,p)

  % convert values to fit into the range [-1,1]

  normvalues     = [];
  normvalues_std = [];

  if p.omit_zeros, nonzero = (values~=0); end

  if p.showsign
    if isempty(valuesmax),   valuesmax  =   max(10^-10,max(abs(values(:)))); end
    if isempty(valuesmin),   valuesmin  = - valuesmax   ; end
  else
    if isempty(valuesmax),   valuesmax  = max(values(:)); end
    if isempty(valuesmin),   valuesmin  = min(values(:)); end
  end

  if valuesmin == valuesmax,  valuesmax = 10^-10 + valuesmin; end
  if p.showsign
    % symmetric scale [-M, M], with M given by valuesmin and valuesmax (default: largest absolute value)
    M = max(10^-10, max(abs([valuesmin, valuesmax])));
    normvalues     = values / M;
    normvalues_std = real(values_std / M);
    valuesmax      = M;
    valuesmin      = -M;
  else
    normvalues     = 2 * [(values - valuesmin)   / ( valuesmax - valuesmin)] - 1;
    normvalues_std = real(values_std /  ( valuesmax- valuesmin));
  end

  normvalues(find(normvalues>1)) = 1;
  normvalues(find(normvalues<-1))=-1;

  if p.omit_zeros
    normvalues = normvalues .* nonzero;
    if ~isempty(normvalues_std), normvalues_std = normvalues_std .* nonzero;  end
  end

  if isnan(valuesmin)
    valuesmin = -abs(valuesmax);
  end

  % ----------------------------------------

function this_colorbar(ax,num,colmap,location,colorbar_fontsize)

  % colorbar outside the axes (location 'South' means 'southoutside', etc.)
  colormap(ax, colmap);
  location = lower(location);
  if isempty(strfind(location,'outside')), location = [location 'outside']; end
  h = colorbar(ax, 'Location', location);
  limits = [min(num), max(num)];
  caxis(ax, limits);
  set(h, 'Limits', limits, 'FontSize', colorbar_fontsize);

  % ----------------------------------------

function set_axis_limits(ax, x, squaresize, keep_old_limits)

  % axis limits that contain all nodes, symbols, and text labels (with a small margin)

  x = x(:, all(isfinite(x),1));
  if isempty(x), return; end
  margin = max(squaresize, 0.02 * max(max(x,[],2) - min(x,[],2)));
  if margin == 0, margin = 1; end
  box = [min(x(1,:))-margin, max(x(1,:))+margin, min(x(2,:))-margin, max(x(2,:))+margin];
  if keep_old_limits
    a = axis(ax);
    box = [min(box(1),a(1)), max(box(2),a(2)), min(box(3),a(3)), max(box(4),a(4))];
  end
  axis(ax, box);

  labels = findobj(ax, 'Type', 'text');
  if isempty(labels), return; end
  [anchor, cx, cy] = label_boxes(labels);
  if isempty(anchor), return; end

  % size of the plot box in points
  units = get(ax, 'Units');
  set(ax, 'Units', 'points');
  plot_box = get(ax, 'Position');
  set(ax, 'Units', units);
  plot_box = max(plot_box(3:4), 1);

  % enlarge the limits to include the labels (iterated, because the size of
  % the labels in data units changes with the limits)
  a = box;
  for pass = 1:3
    scale = max((a(2)-a(1)) / plot_box(1), (a(4)-a(3)) / plot_box(2));   % data units per point
    ex = anchor(1,:) + scale * cx;
    ey = anchor(2,:) + scale * cy;
    a = [min([box(1), ex(:)']), max([box(2), ex(:)']), min([box(3), ey(:)']), max([box(4), ey(:)'])];
    a = a + 0.01 * max(a(2)-a(1), a(4)-a(3)) * [-1 1 -1 1];
  end
  axis(ax, a);


  % -------------------------------------------------

function [anchor, cx, cy] = label_boxes(labels)

  % positions of text labels (2 x n, data units) and corners of the label boxes relative to
  % these positions (4 x n, in points). The size of a label is estimated from its font size and
  % number of characters, which is much faster than asking Matlab for the 'Extent' of each label

  props = get(labels(:), {'Position','String','FontSize','Rotation','HorizontalAlignment','VerticalAlignment','Units'});
  keep  = strcmp(props(:,7), 'data') & ~cellfun(@isempty, props(:,2));
  props = props(keep,:);
  n = size(props,1);
  anchor = zeros(2,n);
  cx = zeros(4,n);
  cy = zeros(4,n);
  for k = 1:n
    [position, str, font_size, rotation, h_align, v_align] = props{k,1:6};
    lines = cellstr(str);
    lines = regexp(strjoin(lines(:)', newline), '\n', 'split');
    lines = regexprep(lines, '\\[a-zA-Z]+', 'x');   % TeX symbols such as \pm count as one character
    lines = regexprep(lines, '[_^{}]', '');
    w = 0.6 * font_size * max(cellfun(@length, lines));
    h = 1.2 * font_size * numel(lines);
    switch h_align
      case 'center', x0 = -w/2;
      case 'right',  x0 = -w;
      otherwise,     x0 = 0;
    end
    switch v_align
      case {'top','cap'},          y0 = -h;
      case {'bottom','baseline'},  y0 = 0;
      otherwise,                   y0 = -h/2;
    end
    corners = [x0, x0+w, x0+w, x0; y0, y0, y0+h, y0+h];
    rot = [cosd(rotation), -sind(rotation); sind(rotation), cosd(rotation)];
    corners = rot * corners;
    anchor(:,k) = position(1:2)';
    cx(:,k) = corners(1,:)';
    cy(:,k) = corners(2,:)';
  end

  % -------------------------------------------------


function p = netgraph_draw_set_parameters(network,vararg)

  [nm,nr] = size(network.N);

  % use parameters from network

  if isfield(network,'graphics_par')
    p = network.graphics_par;
  else
    network = netgraph_make_graph(network);
    p = network.graphics_par;
  end

  % -------------------------
  % use parameters in argument list

  input_arg = netgraph_options(vararg, network);
  fn = fieldnames(input_arg);
  for i = 1:length(fn)
    p.(fn{i}) = input_arg.(fn{i});
  end

  % -------------------------
  % add missing parameters

  p_default = struct(...
      'N',               network.N, ...
      'show_regulation', 0, ...
      'showsign',        1, ...
      'omit_zeros',      0, ...
      'show_metvalues', 1, ...
      'show_actvalues', 1, ...
      'show_arrowvalues', 1, ...
      'show_edgevalues', 0, ...
      'show_regulationvalues', 0, ...
      'metvalues_std',   zeros(size(p.metvalues)),  ...
      'metvaluesmax',    [], ...
      'metvaluesmin',    [], ...
      'metstyle', 'fixed', ...
      'metcolors', [], ...
      'metinvisible', zeros(size(network.N,1),1), ...
      'circle_shift',    [0,0], ...
      'actvalues_std',   zeros(size(p.actvalues)),  ...
      'actvaluesmax',    [], ...
      'actvaluesmin',    [], ...
      'actstyle', 'fixed', ...
      'actcolors', [], ...
      'actinvisible', zeros(size(network.N,2),1), ...
      'omitreactions',   {'omitthisreaction'}, ...
      'omitmetabolites', {'omitthismetabolite'}, ...
      'suppress_lines', 0, ...
      'arrowstyle',      'none', ...
      'arrowvalues',     [], ...
      'arrowprintvalues',   0, ...
      'arrowvaluesmax',  [], ...
      'arrowvaluesmin',  [], ...
      'arrow_stoichiometries', 0, ...
      'arrowvalues_std',     zeros(size(p.arrowvalues)), ...
      'arrow_shift',     0.7,  ...
      'arrowcolor',      [0.55 0.55 0.55],  ...
      'colorbar_fontsize', 10, ...
      'linewidth',       1, ...
      'single_arrow',    0, ...
      'edgevalues',      [],  ...
      'edgevalues_std',  [], ...
      'edgevaluesmax',   [], ...
      'edgevaluesmin',   [], ...
      'edgeprintvalues',      [],  ...
      'edgestyle',      'normal',  ...
      'regulationvalues',    [],  ...
      'regulation_matrix',    [],  ...
      'regulationprintvalues',    [],  ...
      'regulationvalues_std',[], ...
      'regulationvaluesmax', [], ...
      'regulationvaluesmin', [], ...
      'regulationstyle',      'normal',  ...
      'hold_on',         0, ...
      'FontSize',        8, ...
      'fontangle',       0, ...
      'YDir',            'normal', ...
      'figure_axis',     [], ...
      'black_and_white', 0, ...
      'colormap',        rb_colors(250), ...
      'colorbar_location', 'South', ...
      'canvas_position', [], ...
      'keep_subplot', 0, ...
      'linecolor',       [0.3 0.3 0.3], ...
      'rlinecolor',      [1 .3 0.], ...
      'text_offset',     [], ...
      'value_text_offset',     [0.01,-0.01], ...
      'straightlines',   1,       ...
      'arrowsize',       0.03,        ...
      'figure_position', [],  ...
      'subplot_position',[], ...
      'no_points',       0, ...
      'colorbar',        0,...
      'colorbar_numbers', [],...
      'flag_edges', 1, ...
      'flag_triangle_edges', 0, ...
      'canvas', [],...
      'scale_sizes', 1, ...
      'shade_long_lines', [], ...
      'shade_cofactor_lines', []);

  p = join_struct(p_default,p);

  % -------------------------

  if strcmp(p.arrowstyle,'none')
    p.single_arrow = 0;
  end

  if p.single_arrow
    p.show_arrowvalues = 1;
    p.arrow_shift = 1;
  end

  if ~isempty(p.edgevalues),       p.show_edgevalues = 1; end
  if ~isempty(p.regulationvalues), p.show_regulationvalues = 1; end
  if p.show_regulationvalues,      p.show_regulation = 1; end

  if p.arrowprintvalues, p.actprintvalues = 0; end

  % ----------
  % correct parameters

  switch p.arrowstyle
    case 'directions'
      p.arrowvalues = ones(length(network.actions),1);
      if isfield(p,'reaction_mapping')
        p.arrowvalues = ones(max(p.reaction_mapping),1);
      end
  end

  if isempty(p.metvalues),          p.show_metvalues = 0; end
  if strcmp(p.metstyle,'none'),     p.show_metvalues = 0; end
  if isempty(p.actvalues),          p.show_actvalues = 0; end
  if strcmp(p.actstyle,'none'),     p.show_actvalues = 0; end
  if isempty(p.arrowvalues),        p.show_arrowvalues = 0; end
  if isempty(p.edgevalues),         p.show_edgevalues = 0; end
  if isempty(p.regulationvalues),   p.show_regulationvalues = 0; end
  if strcmp(p.arrowstyle,'none'),   p.show_arrowvalues = 0; end

  if p.show_regulationvalues
    if isempty(p.regulationvalues_std)
      p.regulationvalues_std = zeros(size(p.regulationvalues));
    end
  end

  if p.show_edgevalues
    if isempty(p.edgevalues_std)
      p.edgevalues_std = zeros(size(p.edgevalues));
    end
  end

  % ----------
  % split metabolites?

  if isfield(p,'metabolite_mapping')
    if p.show_metvalues
      p.metvalues       = p.metvalues(p.metabolite_mapping);
      p.metvalues_std   = p.metvalues_std(p.metabolite_mapping);
    end
    if p.show_edgevalues
      p.edgevalues     = p.edgevalues(p.metabolite_mapping,:);
      p.edgevalues_std = p.edgevalues_std(p.metabolite_mapping,:);
    end
    if p.show_regulationvalues
      p.regulationvalues     = p.regulationvalues(:,p.metabolite_mapping);
      p.regulationvalues_std = p.regulationvalues_std(:,p.metabolite_mapping);
    end
  end

  if isfield(p,'reaction_mapping')
    if p.show_actvalues
      p.actvalues       = p.actvalues(p.reaction_mapping);
      p.actvalues_std   = p.actvalues_std(p.reaction_mapping);
    end
    if p.show_arrowvalues
      p.arrowvalues       = p.arrowvalues(p.reaction_mapping);
    end
    if p.show_edgevalues
      p.edgevalues     = p.edgevalues(:,p.reaction_mapping);
      p.edgevalues_std = p.edgevalues_std(:,p.reaction_mapping);
    end
    if p.show_regulationvalues
      p.regulationvalues     = p.regulationvalues(p.reaction_mapping,:);
      p.regulationvalues_std = p.regulationvalues_std(p.reaction_mapping,:);
    end
  end

  % ---------------------------------------------------------------------
  % threshold the values at given max and min

  if p.colorbar
    if length(p.colorbar_numbers)
      p.metvaluesmin  = min(p.colorbar_numbers);
      p.metvaluesmax  = max(p.colorbar_numbers);
      p.actvaluesmin  = min(p.colorbar_numbers);
      p.actvaluesmax  = max(p.colorbar_numbers);
      p.edgevaluesmin = min(p.colorbar_numbers(:));
      p.edgevaluesmax = max(p.colorbar_numbers(:));
    end
  end

  if isempty(p.metvaluesmin), p.metvaluesmin = min(p.metvalues); end
  if isempty(p.metvaluesmax), p.metvaluesmax = max(p.metvalues); end
  if isempty(p.actvaluesmin), p.actvaluesmin = min(p.actvalues); end
  if isempty(p.actvaluesmax), p.actvaluesmax = max(p.actvalues); end
  if isempty(p.edgevaluesmin), p.edgevaluesmin = min(p.edgevalues(:)); end
  if isempty(p.edgevaluesmax), p.edgevaluesmax = max(p.edgevalues(:)); end

  p.metvalues( find(p.metvalues < p.metvaluesmin) ) = p.metvaluesmin;
  p.metvalues( find(p.metvalues > p.metvaluesmax) ) = p.metvaluesmax;
  p.actvalues( find(p.actvalues < p.actvaluesmin) ) = p.actvaluesmin;
  p.actvalues( find(p.actvalues > p.actvaluesmax) ) = p.actvaluesmax;

  if size(p.edgevalues)
    p.edgevalues(find(p.edgevalues < p.edgevaluesmin)) = p.edgevaluesmin;
    p.edgevalues(find(p.edgevalues > p.edgevaluesmax)) = p.edgevaluesmax;
  end

  % ---------------------------------------------------------------------
  % showsign = 1  -> normalise all values to values between -1 and 1
  % showsign = 0  -> normalise all values to values between  0 and 1

  if p.show_metvalues
    [p.norm_metvalues, p.norm_metvalues_std,p.metvaluesmax,p.metvaluesmin] = normalise(p.metvalues,p.metvalues_std, p.metvaluesmax, p.metvaluesmin, p);
  end

  if p.show_actvalues
    [p.norm_actvalues, p.norm_actvalues_std,p.actvaluesmax,p.actvaluesmin] = normalise(p.actvalues,p.actvalues_std,p.actvaluesmax,p.actvaluesmin,p);
  end

  if p.show_arrowvalues
    p.norm_arrowvalues = p.arrowvalues;
    if length(p.arrowvaluesmax)
      if isempty(p.arrowvaluesmin)
        p.arrowvaluesmin = -p.arrowvaluesmax;
      end
      p.norm_arrowvalues(p.norm_arrowvalues>p.arrowvaluesmax) = p.arrowvaluesmax;
      p.norm_arrowvalues(p.norm_arrowvalues<p.arrowvaluesmin) = p.arrowvaluesmin;
    end
  else
    p.norm_arrowvalues = zeros(nr,1);
  end

  if p.show_edgevalues
    [p.norm_edgevalues, p.norm_edgevalues_std,p.edgevaluesmax,p.edgevaluesmin] = normalise(p.edgevalues,p.edgevalues_std,p.edgevaluesmax,p.edgevaluesmin,p);
  end

  if p.show_regulationvalues
    [p.norm_regulationvalues, p.norm_regulationvalues_std,p.regulationvaluesmax,p.regulationvaluesmin] = normalise(p.regulationvalues,p.regulationvalues_std,p.regulationvaluesmax,p.regulationvaluesmin,p);
  end

  % --- some more stuff

  if p.black_and_white,        p.colormap = 0.5-0.5*sign(gray(250)-0.5); end

  % omit metabolites that are not connected anymore once reactions have been omitted

  ll = setdiff(1:nr,label_names(p.omitreactions,network.actions));
  p.omitmetabolites = union(p.omitmetabolites, network.metabolites(find(sum(abs(network.N(:,ll)),2)==0)));
  p.actindshow = setdiff(1:length(network.actions),    label_names(p.omitreactions,network.actions));
  p.metindshow = setdiff(1:length(network.metabolites),label_names(p.omitmetabolites,network.metabolites));

  % omit metabolites and reactions that are tagged as invisible
  p.actindshow = p.actindshow(find(p.actinvisible(p.actindshow)==0));
  p.metindshow = p.metindshow(find(p.metinvisible(p.metindshow)==0));

  % ----------

  if isfield(p,'chop_names')
    for it = 1:length(p.metnames)
      if length(p.metnames{it}) > p.chop_names
        p.metnames{it} = [p.metnames{it}(1:p.chop_names) '...'];
      end
    end
    for it = 1:length(p.actnames)
      if length(p.actnames{it}) > p.chop_names
        p.actnames{it} = [p.actnames{it}(1:p.chop_names) '...' ];
      end
    end
  end

  if isfield(p,'split_names')
    if p.split_names
      for it = 1:length(p.metnames)
        p.metnames{it} =  split_name(p.metnames{it},p.split_names);
      end
      for it = 1:length(p.actnames)
        p.actnames{it} =  split_name(p.actnames{it},p.split_names);
      end
    end
  end

  if size(p.x,2) == 0, error('mnt:netgraph:noPositions', 'No node positions found'); end

  % text labels to the right of the symbols (offsets for a layout of size 1)
  if isempty(p.text_offset)
    p.text_offset = [0.5 * p.squaresize + 0.008, 0];
  end

  % size of the layout (largest extent of the shown nodes); symbol sizes are scaled by this number
  p.layout_size = 1;
  if p.scale_sizes
    x_shown = p.x(:, [p.metindshow, nm + p.actindshow]);
    x_shown = x_shown(:, all(isfinite(x_shown),1));
    if size(x_shown,2) > 1
      extent = max(max(x_shown,[],2) - min(x_shown,[],2));
      if extent > 0, p.layout_size = extent; end
    end
  end
