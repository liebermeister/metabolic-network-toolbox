function M_time = netgraph_flux_arrows_movie(network, J_matrix, text_flag, varargin)

% NETGRAPH_FLUX_ARROWS_MOVIE Movie of flux distributions shown by arrows
%
% M_time = netgraph_flux_arrows_movie(network, J_matrix, text_flag, options)
% M_time = netgraph_flux_arrows_movie(network, J_matrix, text_flag, 'name1', value1, ...)
% M_time = netgraph_flux_arrows_movie(network, J_matrix, options)          (text_flag = 0)
%
% J_matrix   matrix of fluxes (one column per frame)
% text_flag  show names of metabolites and reactions (default 0)
% options    graphics options (struct or name-value pairs; see 'help netgraph_draw'), plus
%              timebar            (default 0)
%              background_colors  (default []) RGB colours, interpolated over the frames
%
% See also: netgraph_flux_movie, netgraph_movie

if ~exist('text_flag','var'), text_flag = 0; end

% options given in place of text_flag
if isstruct(text_flag) || ischar(text_flag) || isstring(text_flag)
  varargin  = [{text_flag}, varargin];
  text_flag = 0;
end

options = netgraph_options(varargin, network);

options = join_struct(struct('timebar',0,'background_colors',[]), options);

J_matrix_scaled = J_matrix / max(abs(J_matrix(:)));

clear M_time
ns          = size(J_matrix,2);
time_points = cellstr(num2str([1:ns]'));

options.background_colors = netgraph_movie_background_colors(options.background_colors, 1, ns);

for it = 1:ns
  options.arrowstyle  = 'fluxes';
  options.arrowvalues = J_matrix_scaled(:,it);
  options.arrow_shift = it/(ns+1);
  options.arrowcolor = [1 0 0];
  subplot('position',[ 0 0 1 1]);
  if ~isempty(options.background_colors), options.canvas = options.background_colors(it,:); end
  netgraph_concentrations(network,[],[],text_flag,options);
  if it==1, a = axis; end
  if options.timebar
    ss = network.graphics_par.squaresize;
    netgraph_movie_timebar(a, ss, (it-1)/max(ns-1,1), ss/2);
    axis(a);
    text(0.01,0.02,time_points{it});
  end
  drawnow
  M_time(it) = getframe;
end
