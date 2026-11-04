% NETGRAPH_MOVIE Create a movie of a network with changing values
%
% [M,T] = netgraph_movie(network, t, s_t, data_type, n_frames, text_flag, gp)
% [M,T] = netgraph_movie(network, t, s_t, data_type, n_frames, text_flag, 'name1', value1, ...)
%
% s_t must be #variables x #timepoints !!!
%
% data_type {'concentrations','reactions','both','compute reactions'}
% default: concentrations
%
% if 'both' is chosen, then s_t must contain the matrix of concentration
% time courses and below the matrix of reaction velocity time courses
%
% Symbol colours for each frame can be specified in gp.metabolite_colors and gp.reaction_colors
%
% usage: movie(M);
%
% gp: graphics gp (struct or name-value pairs; see 'help netgraph_draw'), plus these fields
% (with default values):
%   flux_moves            0
%   metabolite_colors     []
%   reaction_colors       []
%   metvaluesmin          0
%   timebar               1
%   background_colors     []
%   use_background_colors 0
%   rest_frames
%   prefix
%
% to save a movie as an animated gif, use movie_save(filename,M)
%
% See also: netgraph_flux_movie, netgraph_draw

function [M,T] = netgraph_movie(network,t,s_t,data_type,n_frames,text_flag,varargin)

eval(default('data_type','''concentrations''','n_frames','[]','text_flag','0'));

gp = netgraph_options(varargin, network);

gp_default = struct('flux_moves',0,'metabolite_colors',[],'reaction_colors',[],'metvaluesmin',0,'timebar',1,'background_colors',[],'use_background_colors',0);
gp_default.timebar_labels = {};

gp = join_struct(gp_default, gp);

if isempty(n_frames), n_frames = 10; end

t = t-min(t);

T = max(t)*(0:1/(n_frames-1):1);

nandata = find(sum(isnan(s_t),2));
s_t(nandata,:)=0;
s_t            = interp1(t,s_t',T,'pchip')';
s_t(nandata,:) = nan;
if size(s_t,2)==1, s_t = s_t'; end

if length(gp.metabolite_colors)
  gp.metabolite_colors = interp1(t,gp.metabolite_colors,T,'pchip');
  gp.metabolite_colors(gp.metabolite_colors<0) = 0;
  gp.metabolite_colors(gp.metabolite_colors>1) = 1;
end

if length(gp.reaction_colors)
  gp.reaction_colors = interp1(t,gp.reaction_colors,T,'pchip');
  gp.reaction_colors(gp.reaction_colors<0) = 0;
  gp.reaction_colors(gp.reaction_colors>1) = 1;
end

switch data_type
  case 'concentrations'
    s_t_c = s_t;
    s_t_v = [];
    nm    = size(s_t,1);
    nr    = 0;
    %  case 'compute reactions',
    %    if isfield(network.graphics_par,'metabolite_mapping'),
    %       if max(network.graphics_par.metabolite_mapping) > nm,
    %          error('Option compute reactions is not supported in this case');
    %       end
    %    end
    %    s_t_c = s_t;
    %    s_t_v = [];
    %    nm    = size(s_t,1);
    %    nr    = 0;
  case 'reactions'
    s_t_c = [];
    s_t_v = s_t;
    nm    = 0;
    nr    = size(s_t,1);
  case 'both'
    %% warning('data type "both" is not fully supported');
    [nm,nr] = size(network.N);
    s_t_c = s_t(1:nm,:);
    s_t_v = s_t(nm+1:end,:);
  otherwise
    error('Unknown function option');
end

metvalues = zeros(nm,length(T));
actvalues = zeros(nr,length(T));

for j = 1:length(T)
  switch data_type
    case 'concentrations'
      metvalues(:,j) = s_t_c(:,j);
      gp.actstyle = 'none';
    case 'reactions'
      actvalues(:,j) =  s_t_v(:,j);
      gp.metstyle = 'none';
    case 'both'
      %% warning('data type "both" is not fully supported');
      metvalues(:,j) = s_t_c(:,j);
      actvalues(:,j) = s_t_v(:,j);
    case 'compute reactions'
      error('data type "compute reactions" is currently not supported');
    otherwise, error('Unknown function option');
  end
end

jj = j;

if length(metvalues), gp_default.metvaluesmax = max(max(metvalues(:)),0.001); end
if length(actvalues), gp_default.actvaluesmax = max(max(actvalues(:)),0.001); end

clf;
subplot('position',[0 0 1 1]);
netgraph_concentrations(network,metvalues(:,1),actvalues(:,1),0,gp_default);

gp_default.figure_axis    = axis;
gp_default.actprintvalues = 0;
gp_default.metprintvalues = 0;
gp_default.arrowstyle     = 'none';

gp = join_struct(gp_default,gp);

gp.background_colors = netgraph_movie_background_colors(gp.background_colors, gp.use_background_colors, length(T));

a = axis;
clf;

for j=1:jj

  set(gcf,'Color',[1 1 1]);

  if gp.flux_moves
    arrow_shift         = mod(gp.flux_moves*j,jj)/jj;
    gp.arrow_shift = arrow_shift;
  end

  if length(gp.metabolite_colors)
    gp.metcolors = squeeze(gp.metabolite_colors(j,:,:));
    if size( gp.metcolors,2) == 1, gp.metcolors = gp.metcolors'; end
  end

  if length(gp.reaction_colors)
    gp.actcolors = squeeze(gp.reaction_colors(j,:,:));
    if size( gp.actcolors,2) == 1, gp.actcolors = gp.actcolors'; end
  end

  clf;
  subplot('position',[0 0 1 1]);
  if ~isempty(gp.background_colors)
    gp.canvas = gp.background_colors(j,:);
  end
  netgraph_concentrations(network,metvalues(:,j),actvalues(:,j),text_flag,gp);
  if gp.timebar
    ss = network.graphics_par.squaresize;
    netgraph_movie_timebar(a, ss, (j-1)/max(jj-1,1));
    if length(gp.timebar_labels)
      text(0.9*a(1)+0.1*a(2), 0.05*a(3)+0.95*a(4), gp.timebar_labels{1+floor(length(t)*j/[jj+1])},'FontSize',24);
    end
  end
  axis(a);

  M(j) = getframe(gca);

  if isfield(gp,'prefix')
    if j ==1,             print(gcf,[gp.prefix '_start.eps'],'-depsc'); end
    if j == ceil(0.8*jj), print(gcf,[gp.prefix '_end.eps'],'-depsc');   end
  end

end

if isfield(gp,'rest_frames')
  M = [ repmat(M(1),1,gp.rest_frames), M, repmat(M(end),1,gp.rest_frames)];
end
