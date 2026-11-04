function netgraph_movie_timebar(a, ss, fraction, radius)

% netgraph_movie_timebar(a, ss, fraction, radius)
%
% Draw a time bar at the bottom of a network movie frame (used by the netgraph_..._movie functions)
%
% a         axis limits [xmin xmax ymin ymax]
% ss        distance from the figure border (usually graphics_par.squaresize)
% fraction  position of the marker on the bar (between 0 and 1)
% radius    radius of the marker (default ss/4)

if ~exist('radius','var'), radius = ss/4; end

fraction = min(max(fraction,0),1);
hold on;
plot([a(1)+ss, a(2)-ss], ss + [a(3), a(3)], '-', 'Color', [0 0 0]);
x_marker = a(1) + ss + fraction * (a(2) - a(1) - 2*ss);
circle(x_marker, ss + a(3), radius, [0 0 0]);
circle(x_marker, ss + a(3), ss/6,   [1 1 1]);
