% M = netgraph_flux_movie(network,S,J,text_flag,goptions);
%
% S and J : vectors of concentrations and fluxes
% the number of frames can be set by 'n_frames' within goptions
%
% to save a movie as an animated gif, use movie_save(filename,M)

function M = netgraph_flux_movie(network,S,J,text_flag,goptions);

if ~exist('text_flag','var'),  text_flag=0; end
if ~exist('goptions','var'), goptions=struct; end

arrow_shift = 0;

%if ~isfield(goptions,'arrowsize'), goptions.arrowsize = 0.03; end
if ~isfield(goptions,'n_frames'),  goptions.n_frames  = 10;   end

goptions.actstyle    ='none';
goptions.arrow_shift = arrow_shift;

netgraph_concentrations(network,S,J,text_flag,goptions); a = axis;
netgraph_concentrations(network,S,J,text_flag,goptions); axis(a); axis off;
M(1) = getframe;

for it = 1:goptions.n_frames,
  arrow_shift            = it/(goptions.n_frames+1);
  goptions.arrow_shift   = arrow_shift;
  netgraph_concentrations(network,S,J,text_flag,goptions); 
  axis(a);  axis off;
  M(it+1) = getframe;
end
