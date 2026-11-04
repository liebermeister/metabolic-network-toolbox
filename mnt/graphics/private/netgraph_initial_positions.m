function x = netgraph_initial_positions(graphics_par, indices_fixed, x_fixed)

% x = netgraph_initial_positions(graphics_par, indices_fixed, x_fixed)
%
% Initial guess of node positions for a network layout (used by netgraph_find_x):
% multidimensional scaling of the graph distances, followed by a few steps in which
% nodes repel each other and are attracted by their neighbours
%
% Inputs:
%   graphics_par   network.graphics_par, with fields db (graph distances), N, and metnames
%   indices_fixed  indices of nodes (metabolites first, then reactions) with fixed positions
%   x_fixed        positions of these nodes (2 x n matrix)
%
% Output:
%   x  node positions (2 x number of nodes), scaled to the range [0,1] and rounded to a grid
%
% (Formerly netgraph_find_x_old in matlab-utils.)

if isempty(graphics_par.db)
  x = [];
  return
end

n_iterations = 12;
gridsize = 0.01;                                     % distance between discrete x values
k1 = 0.005;                                          % overall repulsion strength
k2 = 0.5 * sqrt(1/length(graphics_par.metnames));    % overall repulsion scale
k3 = 0.2;                                            % neighbour attraction strength

db = graphics_par.db.^0.2;
x  = multidimensional_scaling(db);
x  = diag([-1 1]' .* sign(x(:,1) - mean(x,2))) * x;

[Ni,Nj] = ind2sub(size(graphics_par.N), find(graphics_par.N ~= 0));
n_nodes = size(x,2);

% move points around -> similar distances
for it2 = 1:n_iterations
  x = x + 0.001 * randn(size(x));
  for it = 1:3
    edgecenters = 0.5 * (x(:,Ni) + x(:,Nj));
    xnew = zeros(size(x));
    for i = 1:n_nodes
      distances = repmat(x(:,i), 1, n_nodes + length(Ni)) - [x, edgecenters];
      f1 = sum(distances.^2) + 10^-10;
      f1 = exp(-f1/(2*k2^2)) ./ sqrt(f1);
      neighbours = x(:, db(i,:) == 1);
      if isempty(neighbours)
        x_neighbours = x(:,i);
      else
        x_neighbours = mean(neighbours,2);
      end
      xnew(:,i) = x(:,i) + k1 * sum(distances .* repmat(f1,2,1), 2) + k3 * (x_neighbours - x(:,i));
    end
    x = xnew;
    if ~isempty(indices_fixed)
      x(:,indices_fixed) = x_fixed;
    end
  end
end

x = x - repmat(min(x,[],2), 1, n_nodes);
x = diag(1 ./ max(x,[],2)) * x;
x = gridsize * round(x/gridsize);
