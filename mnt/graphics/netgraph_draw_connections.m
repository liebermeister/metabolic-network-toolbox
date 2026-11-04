% NETGRAPH_DRAW_CONNECTIONS Draw lines between reactions, for a matrix of relations between reactions
%
% netgraph_draw_connections(network, M)
%
% M  symmetric matrix (reactions x reactions); for non-zero values, a line is drawn
%    between the two reactions (blue for negative values, red for positive values)
%
% The lines are added to the current axes (e.g. after netgraph_draw).
%
% See also: netgraph_draw, synergy_network_plot

function netgraph_draw_connections(network,M)

n_react = length(network.graphics_par.actnames);
x=network.graphics_par.x(:,end-n_react+1:end);
for i1 = 1: n_react
  for i2 = 1: i1-1
    % blue for negative, red for positive values
    c = [0 0 1]*(M(i1,i2)<0) + [1 0 0]*(M(i1,i2)>0);
    if M(i1,i2)~=0
      h=line( [x(1,i1) x(1,i2)],[x(2,i1) x(2,i2)]);
      set(h,'Color',c);
    end
  end
end
