function M_time = netgraph_flux_arrows_movie(network, J_matrix, text_flag, gp)

% NETGRAPH_FLUX_ARROWS_MOVIE Movie of flux distributions shown by arrows
%
% M_time = netgraph_flux_arrows_movie(network, J_matrix, text_flag, gp)
%
% J_matrix   matrix of fluxes (one column per frame)
% text_flag  show element names (default 0)
% gp         graphics options (see netgraph_draw), with extra fields
%              timebar            (default 0)
%              background_colors  (default []) RGB colours, interpolated over the frames
%
% See also: netgraph_flux_movie, netgraph_movie

eval(default('text_flag','0','gp','struct'));

gp = join_struct(struct('timebar',0,'background_colors',[]), gp);

J_matrix_scaled = J_matrix / max(abs(J_matrix(:)));

clear M_time
ns          = size(J_matrix,2);
time_points = cellstr(num2str([1:ns]'));

gp.background_colors = netgraph_movie_background_colors(gp.background_colors, 1, ns);

for it = 1:ns
  gp.arrowstyle  = 'fluxes';
  gp.arrowvalues = J_matrix_scaled(:,it);
  gp.arrow_shift = it/(ns+1);
  gp.arrowcolor = [1 0 0];
  subplot('position',[ 0 0 1 1]);
  if ~isempty(gp.background_colors), gp.canvas = gp.background_colors(it,:); end
  netgraph_concentrations(network,[],[],text_flag,gp);
  if it==1, a = axis; end
  if gp.timebar
    ss = network.graphics_par.squaresize;
    netgraph_movie_timebar(a, ss, (it-1)/max(ns-1,1), ss/2);
    axis(a);
    text(0.01,0.02,time_points{it});
  end
  drawnow
  M_time(it) = getframe;
end
