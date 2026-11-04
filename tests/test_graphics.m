% Tests for network graphics (folder 'graphics')

setup_test_path;
network = netgraph_make_graph(example_network);
nm = 6; nr = 5;

%% netgraph_draw: default settings
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
netgraph_draw(network);
assert(~isempty(findobj(gca,'Type','line')));

%% netgraph_draw: curved lines, with and without arrows
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
netgraph_draw(network, struct('straightlines',0,'arrowstyle','directions'));
assert(size(tagged_xdata('netgraph_arrows'),2) == nnz(network.N));
clf; netgraph_draw(network, struct('straightlines',0,'arrowstyle','none'));
% one line per edge (lines are separated by nan values), no arrows
assert(sum(isnan(tagged_xdata('netgraph_edges'))) == nnz(network.N));
assert(isempty(tagged_xdata('netgraph_arrows')));

%% netgraph_draw: metvaluesmax sets the scale of the symbols
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
gp = struct('metvalues',(1:6)'/10,'metstyle','box','metvaluesmin',0,'arrowstyle','none', ...
            'actstyle','none','metprintnames',0,'actprintnames',0);
width = zeros(1,2); maxvalues = [0.6, 1.2];
for it = 1:2
  clf; netgraph_draw(network, join_struct(gp, struct('metvaluesmax',maxvalues(it))));
  width(it) = max(symbol_widths('netgraph_metabolites'));
end
assert(abs(width(2)/width(1) - 0.5) < 0.05);

%% netgraph_draw: with showsign = 1, the colorbar is symmetric
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
netgraph_draw(network, struct('metvalues',(1:6)','colorbar',1,'arrowstyle','none'));
cb = findobj(fig,'Type','ColorBar');
assert(abs(cb.Limits(1) + cb.Limits(2)) < 10^-10 && cb.Limits(2) >= 6);

%% netgraph_draw: single arrows for positive and negative fluxes
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
for flux_sign = [1, -1]
  clf; netgraph_draw(network, struct('arrowstyle','fluxes','arrowvalues',flux_sign*ones(nr,1), ...
                                     'single_arrow',1,'actstyle','none','metprintnames',0,'actprintnames',0));
  assert(size(tagged_xdata('netgraph_arrows'),2) == nr);
end

%% netgraph_move: display settings are not stored in the network
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
moved = netgraph_move(network, 'none', struct('manual',0));
assert(strcmp(moved.graphics_par.arrowstyle, network.graphics_par.arrowstyle));

%% netgraph_make_graph: the random number generator is not changed
rng(5); a = rand;
rng(5); netgraph_make_graph(example_network); b = rand;
assert(a == b);

%% netgraph_make_graph: fixed positions in networks with more than 255 nodes
N = zeros(150,149);
for k = 1:149, N(k,k) = -1; N(k+1,k) = 1; end
n_nodes = 150 + 149;
x_fixed = round(100*rand(2,n_nodes))/100;
x_fixed(:,1) = [0;0]; x_fixed(:,2) = [1;1];
big = netgraph_make_graph(network_construct(N), [], [], (1:n_nodes)', x_fixed);
% positions are only rounded to a fine grid (spacing 0.01/sqrt(150))
assert(max(abs(big.graphics_par.x(:)-x_fixed(:))) < 10^-3);

%% netgraph_print_positions and netgraph_read_positions: elements missing in model or file
filename = [tempname '.tsv'];
cleanup = onCleanup(@() delete(filename));
netgraph_print_positions(network, filename, [0 0], 'replace file', 0);
other = network;
other.metabolites{1} = 'new_metabolite';
other.graphics_par.fixed(:) = 0; other.graphics_par.fixed(2) = 1;
[names, positions] = netgraph_print_positions(other, filename, [0 0], 'replace elements', 0);
t = sbtab_table_load(filename);
elements = sbtab_table_get_column(t,'Element');
assert(numel(elements) == nm + nr + 1);
assert(all(ismember([other.metabolites; other.actions], elements)));
netgraph_print_positions(other, filename, [0 0], 'add nonexistent', 0);
again = netgraph_read_positions(other, filename, [0 0], 1, 0);
assert(max(abs(again.graphics_par.x(:) - other.graphics_par.x(:))) < 10^-6);

%% netgraph_split_metabolites: network without regulation matrix
plain = network;
if isfield(plain,'regulation_matrix'), plain = rmfield(plain,'regulation_matrix'); end
split = netgraph_split_metabolites(plain, {'S3'});
assert(numel(split.metabolites) == nm + 2);

%% netgraph_flux_movie and netgraph_movie: time bar and background colours
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
lastwarn('');
M = netgraph_flux_movie(network, [], ones(nr,1), 0, ...
                        struct('timebar',1,'n_frames',3,'use_background_colors',1,'background_colors',[1 0 0; 0 0 1]));
assert(numel(M) == 4);
assert(isempty(lastwarn));
M = netgraph_movie(network, [0 1 2], rand(nm,3), 'concentrations', 4, 0, ...
                   struct('use_background_colors',1,'background_colors',[1 1 1; 1 0 0]));
assert(numel(M) == 4);

%% synergy_network_plot: values for pairs of reactions
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
pair_values = [0 1 0 0 0; 1 0 -1 0 0; 0 -1 0 2 0; 0 0 2 0 0; 0 0 0 0 0];
synergy_network_plot(network, ones(nr,1), pair_values);
assert(numel(findobj(gca,'Type','line')) >= 3);

%% netgraph_concentrations in subplots: both subplots are kept
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
subplot(1,2,1); netgraph_concentrations(network, rand(nm,1), [], 1);
subplot(1,2,2); netgraph_concentrations(network, [], rand(nr,1), 1);
assert(numel(findobj(fig,'Type','axes')) == 2);

%% netgraph_draw: symbol sizes scale with the size of the layout
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
gp = struct('metvalues',ones(nm,1),'arrowstyle','none','actstyle','none','metprintnames',0,'actprintnames',0);
netgraph_draw(network, gp);
w1 = max(symbol_widths('netgraph_metabolites'));
big = network; big.graphics_par.x = 500 * big.graphics_par.x;
clf; netgraph_draw(big, gp);
w2 = max(symbol_widths('netgraph_metabolites'));
assert(abs(w2/w1 - 500) < 1);

%% netgraph_draw: all labels lie within the axis limits
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
netgraph_draw(network);
a = axis;
e = get(findobj(gca,'Type','text'), {'Extent'}); e = vertcat(e{:});
assert(all(e(:,1) >= a(1)) && all(e(:,1)+e(:,3) <= a(2)) && all(e(:,2) >= a(3)) && all(e(:,2)+e(:,4) <= a(4)));

%% netgraph_read_positions: elements missing in the file are placed near their neighbours
filename = [tempname '.tsv'];
cleanup = onCleanup(@() delete(filename));
netgraph_print_positions(network, filename, [0 0], 'replace file', 0);
other = network;
other.metabolites{3} = 'new_metabolite';
other.actions{3}     = 'new_reaction';
other = netgraph_read_positions(other, filename, [0 0], 1, 0);
x = other.graphics_par.x;
assert(all(isfinite(x(:))));
assert(norm(x(:,3) - x(:,nm+3)) > 0);
% the new reaction lies near its metabolite S4
assert(norm(x(:,nm+3) - x(:,4)) < norm(network.graphics_par.x(:,1) - network.graphics_par.x(:,6)));

%% netgraph_make_graph: names shown in graphics
n = network_construct([-1; 1]);
n.metabolites = {'Acetyl-CoA'; 'UDP-D-glucose'};
n = netgraph_make_graph(n);
assert(isequal(n.graphics_par.metnames, {'Acetyl-CoA'; 'UDP-glucose'}));

%% netgraph_draw: unknown fluxes (nan) are shown by grey diamonds
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
v = [1 nan 1 1 1]';
for single_arrow = [0 1]
  for straightlines = [1 0]
    clf; netgraph_draw(network, struct('arrowstyle','fluxes','arrowvalues',v,'single_arrow',single_arrow, ...
                                       'straightlines',straightlines,'actstyle','none'));
    % one diamond per edge of reaction 2 (S2-E2, E2-S3), or one diamond with single_arrow
    assert(size(tagged_xdata('netgraph_unknown_fluxes'),2) == (single_arrow == 0) * nnz(network.N(:,2)) + (single_arrow == 1));
    % arrows for the other reactions
    assert(size(tagged_xdata('netgraph_arrows'),2) == (single_arrow == 0) * (nnz(network.N) - nnz(network.N(:,2))) + (single_arrow == 1) * (nr-1));
  end
end

%% netgraph_make_graph: positions from a layout file
filename = [tempname '.tsv'];
cleanup = onCleanup(@() delete(filename));
netgraph_print_positions(network, filename, [0 0], 'replace file', 0);
n = netgraph_make_graph(example_network, [], [], [], [], filename);
assert(max(abs(n.graphics_par.x(:) - network.graphics_par.x(:))) < 10^-6);

%% Graphics options: struct, name-value pairs, and options in place of text_flag give the same result
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
c = (1:nm)'; v = ones(nr,1);
netgraph_concentrations(network, c, v, 1, struct('metprintnames',0,'actprintnames',1,'FontSize',12));
n1 = numel(findobj(gca,'Type','text'));
clf; netgraph_concentrations(network, c, v, 1, 'metprintnames', 0, 'actprintnames', 1, 'FontSize', 12);
n2 = numel(findobj(gca,'Type','text'));
clf; netgraph_concentrations(network, c, v, struct('colorbar',1));
assert(n1 == n2 && n1 == nr);
assert(~isempty(findobj(fig,'Type','ColorBar')));
assert(isempty(findobj(gca,'Type','text')));
clf; netgraph_concentrations(network, c, v, 'colorbar', 1);
assert(~isempty(findobj(fig,'Type','ColorBar')));

%% Graphics options: upper/lower case is corrected, 'Rotation' is accepted, unknown names give a warning
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
netgraph_draw(network, 'fontsize', 13, 'Rotation', 30, 'metprintnames', 1, 'actprintnames', 0);
t = findobj(gca,'Type','text');
assert(all(cell2mat(get(t,'FontSize')) == 13) && all(cell2mat(get(t,'Rotation')) == 30));
lastwarn('');
clf; netgraph_draw(network, 'no_such_option_for_testing', 1);
[~, id] = lastwarn;
assert(strcmp(id, 'mnt:netgraph:unknownOption'));

%% netgraph_draw: values, standard deviations, edge and regulation symbols
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
netgraph_draw(network, struct('metvalues',[1 2 nan 4 5 6]','metvalues_std',0.1*ones(nm,1),'metstyle','box_std', ...
                              'actvalues',(1:nr)','actstyle','box','actprintvalues',1, ...
                              'edgevalues',sparse(double(network.N~=0)),'arrowstyle','none'));
assert(size(tagged_xdata('netgraph_metabolites'),2) == nm);
assert(size(tagged_xdata('netgraph_metabolites_std'),2) == nm);
assert(size(tagged_xdata('netgraph_reactions'),2) == nr);
assert(size(tagged_xdata('netgraph_edge_symbols'),2) == nnz(network.N));
assert(numel(findobj(gca,'Type','text')) >= nr);
reg = network;
reg.regulation_matrix = sparse(nr,nm); reg.regulation_matrix(1,3) = 1; reg.regulation_matrix(2,5) = -1;
clf; netgraph_draw(reg, 'show_regulation', 1, 'regulationvalues', reg.regulation_matrix);
assert(~isempty(tagged_xdata('netgraph_activation')) && ~isempty(tagged_xdata('netgraph_inhibition')));
assert(size(tagged_xdata('netgraph_regulation_symbols'),2) == 2);

% ------------------------------------------------------------------
% Graphics elements of one kind are drawn as one object (with a tag);
% each column of XData is one symbol

function X = tagged_xdata(tag)
  h = findobj(gca,'Tag',tag);
  X = [];
  if ~isempty(h), X = get(h,'XData'); end
end

function w = symbol_widths(tag)
  X = tagged_xdata(tag);
  w = max(X,[],1) - min(X,[],1);
end
