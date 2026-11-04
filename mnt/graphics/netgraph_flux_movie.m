% NETGRAPH_FLUX_MOVIE Movie of a single flux distribution (moving arrows)
%
% M = netgraph_flux_movie(network, S, J, text_flag, gp)
% M = netgraph_flux_movie(network, S, J, text_flag, 'name1', value1, ...)
% M = netgraph_flux_movie(network, S, J, gp)                  (text_flag = 0)
%
% S          metabolite levels ([] for none)
% J          fluxes
% text_flag  show names of metabolites and reactions (default 0)
% gp         graphics options (struct or name-value pairs; see 'help netgraph_draw'), plus
%              n_frames               (default 10) number of frames
%              timebar                (default 0)
%              background_colors      (default []) RGB colours, interpolated over the frames
%              use_background_colors  (default 0)
%
% To save a movie as an animated gif, use movie_save(filename,M)
%
% See also: netgraph_flux_arrows_movie, netgraph_movie

function M = netgraph_flux_movie(network,S,J,text_flag,varargin)

if ~exist('text_flag','var'), text_flag = 0; end

% options given in place of text_flag
if isstruct(text_flag) || ischar(text_flag) || isstring(text_flag)
  varargin  = [{text_flag}, varargin];
  text_flag = 0;
end

gp = netgraph_options(varargin, network);

% scale flux data
J = 1/max(abs(J(:))) * J;
if isfield(gp,'arrowvalues')
  gp.arrowvalues = 1/max(abs(gp.arrowvalues(:))) * gp.arrowvalues;
else
  gp.arrowvalues = J;
end

arrow_shift = 0;

gp_default = struct('n_frames', 10,'timebar',0,'background_colors',[],'use_background_colors',0,'scale_arrowvalues',0);

gp                = join_struct(gp_default,gp);
gp.actstyle       ='none';
gp.arrow_shift    = arrow_shift;
gp.show_metvalues = 0;
% otherwise arrows may be below lines (matlab error)

gp.background_colors = netgraph_movie_background_colors(gp.background_colors, gp.use_background_colors, gp.n_frames);

subplot('Position', [0 0 1 1]);
netgraph_concentrations(network,S,J,text_flag,gp); a = axis;
subplot('Position', [0 0 1 1]);
netgraph_concentrations(network,S,J,text_flag,gp); axis(a); axis off;
M(1) = getframe;

for it = 1:gp.n_frames
  clf;
  subplot('Position', [0 0 1 1]);
  arrow_shift          = it/(gp.n_frames+1);
  gp.arrow_shift = arrow_shift;
  if ~isempty(gp.background_colors)
    gp.canvas = gp.background_colors(it,:);
  end
  netgraph_concentrations(network,S,J,text_flag,gp);
  if gp.timebar
    ss = network.graphics_par.squaresize;
    netgraph_movie_timebar(a, ss, (it-1)/max(gp.n_frames-1,1), ss/2);
  end
  axis(a); axis off;
  M(it+1) = getframe;
end
