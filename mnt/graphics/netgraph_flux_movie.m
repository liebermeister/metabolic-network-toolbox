% M = netgraph_flux_movie(network,S,J,text_flag,goptions);
%
% Movie showing a single flux distribution (moving arrows only to highlight fluxes)
%
% S and J : vectors of concentrations and fluxes
% the number of frames can be set by 'n_frames' within goptions
%
% to save a movie as an animated gif, use movie_save(filename,M)

function M = netgraph_flux_movie(network,S,J,text_flag,goptions);

eval(default('text_flag','0','goptions','struct'));

% scale flux data
J = 1/max(abs(J(:))) * J;
if isfield(goptions,'arrowvalues')
  goptions.arrowvalues = 1/max(abs(goptions.arrowvalues(:))) * goptions.arrowvalues;
else
  goptions.arrowvalues = J;
end

arrow_shift = 0;

goptions_default = struct('n_frames', 10,'timebar',0,'background_colors',[],'use_background_colors',0,'scale_arrowvalues',0);

goptions                = join_struct(goptions_default,goptions);
goptions.actstyle       ='none';
goptions.arrow_shift    = arrow_shift;
goptions.show_metvalues = 0;
% otherwise arrows may be below lines (matlab error)

goptions.background_colors = netgraph_movie_background_colors(goptions.background_colors, goptions.use_background_colors, goptions.n_frames);

subplot('Position', [0 0 1 1]);
netgraph_concentrations(network,S,J,text_flag,goptions); a = axis;
subplot('Position', [0 0 1 1]);
netgraph_concentrations(network,S,J,text_flag,goptions); axis(a); axis off;
M(1) = getframe;

for it = 1:goptions.n_frames
  clf;
  subplot('Position', [0 0 1 1]);
  arrow_shift          = it/(goptions.n_frames+1);
  goptions.arrow_shift = arrow_shift;
  if ~isempty(goptions.background_colors)
    goptions.canvas = goptions.background_colors(it,:);
  end
  netgraph_concentrations(network,S,J,text_flag,goptions);
  if goptions.timebar
    ss = network.graphics_par.squaresize;
    netgraph_movie_timebar(a, ss, (it-1)/max(goptions.n_frames-1,1), ss/2);
  end
  axis(a); axis off;
  M(it+1) = getframe;
end
