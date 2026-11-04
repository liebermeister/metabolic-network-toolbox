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
assert(numel(findobj(gca,'Type','patch')) == nnz(network.N));
clf; netgraph_draw(network, struct('straightlines',0,'arrowstyle','none'));
% one line per edge, plus the dots at the nodes
assert(numel(findobj(gca,'Type','line')) == nnz(network.N) + 1);

%% netgraph_draw: metvaluesmax sets the scale of the symbols
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
gp = struct('metvalues',(1:6)'/10,'metstyle','box','metvaluesmin',0,'arrowstyle','none', ...
            'actstyle','none','metprintnames',0,'actprintnames',0);
width = zeros(1,2); maxvalues = [0.6, 1.2];
for it = 1:2
  clf; netgraph_draw(network, join_struct(gp, struct('metvaluesmax',maxvalues(it))));
  w = cellfun(@(z) max(z)-min(z), get(findobj(gca,'Type','patch'),'XData'));
  width(it) = max(w);
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
  assert(numel(findobj(gca,'Type','patch')) == nr);
end

%% netgraph_move: display settings are not stored in the network
fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
moved = netgraph_move(network, 'none', struct('manual',0));
assert(strcmp(moved.graphics_par.arrowstyle, network.graphics_par.arrowstyle));

%% netgraph_make_graph: the random number generator is not changed
rng(5); a = rand;
rng(5); netgraph_make_graph(example_network); b = rand;
assert(a == b);

%% netgraph_find_x: fixed positions in networks with more than 255 nodes
N = zeros(150,149);
for k = 1:149, N(k,k) = -1; N(k+1,k) = 1; end
big = netgraph_make_graph(network_construct(N));
n_nodes = 150 + 149;
x_fixed = round(100*rand(2,n_nodes))/100;
x_fixed(:,1) = [0;0]; x_fixed(:,2) = [1;1];
x = netgraph_find_x(big, 0.01, 0.01, 0.1, 0.1, 0.5, [(1:n_nodes)', x_fixed']);
assert(max(abs(x(:)-x_fixed(:))) < 10^-6);

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
