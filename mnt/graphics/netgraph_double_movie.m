% NETGRAPH_DOUBLE_MOVIE Movie of two data sets on a network
%
% [M,T] = netgraph_double_movie(network, t, s_t_1, s_t_2, data_type, n_frames, text_flag, gp)
% [M,T] = netgraph_double_movie(network, t, s_t_1, s_t_2, data_type, n_frames, text_flag, 'name1', value1, ...)
%
% Movie displaying data on network; two data sets are shown at the same time
% (for instance, to compare experimental and simulated time series)
%
% gp: graphics gp (struct or name-value pairs; see 'help netgraph_draw' and 'help netgraph_movie'),
% plus these fields (with default values):
% flux_moves         - 0
% metabolite_colors  - []
% reaction_colors    - []
% set_circle_shift   - [0.01    - 0.01]
% colormap_1         - []
% colormap_2         - []

function [M,T] = netgraph_double_movie(network,t,s_t_1,s_t_2,data_type,n_frames,text_flag,varargin)

eval(default('data_type','''concentrations''','n_frames','[]','text_flag','0'));

gp = netgraph_options(varargin, network);

gp_default = struct('flux_moves',0,'metabolite_colors',[],'reaction_colors',[],'metvaluesmin',0,'timebar',1,'background_colors',[],'use_background_colors',0,'colormap_1',[],'colormap_2',[],'set_circle_shift',[0.01,0.01]);

gp = join_struct(gp_default, gp);
gp.canvas=[];

if isempty(n_frames), n_frames = 10; end

t = t-min(t);

T = max(t)*(0:1/(n_frames-1):1);

nandata = find(sum(isnan(s_t_1),2));
s_t_1(nandata,:)=0;
s_t_1            = interp1(t,s_t_1',T,'pchip')';
s_t_1(nandata,:) = nan;
if size(s_t_1,2)==1, s_t_1 = s_t_1'; end

nandata = find(sum(isfinite(s_t_2),2)==0);
s_t_2(sum(isfinite(s_t_2),2)==0,:)=0;
s_t_2            = interp1(t,s_t_2',T,'pchip')';
s_t_2(nandata,:) = nan;
if size(s_t_2,2)==1, s_t_2 = s_t_2'; end

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
    s_t_1_c = s_t_1;
    s_t_1_v = [];
    s_t_2_c = s_t_2;
    s_t_2_v = [];
    nm    = size(s_t_1,1);
    nr    = 0;
  case 'reactions'
    s_t_1_c = [];
    s_t_1_v = s_t_1;
    s_t_2_c = [];
    s_t_2_v = s_t_2;
    nm    = 0;
    nr    = size(s_t_1,1);
end

metvalues_1 = zeros(nm,length(T));
actvalues_1 = zeros(nr,length(T));
metvalues_2 = zeros(nm,length(T));
actvalues_2 = zeros(nr,length(T));

gp.background_colors = netgraph_movie_background_colors(gp.background_colors, gp.use_background_colors, length(T));


for j = 1:length(T)
  switch data_type
    case 'concentrations'
      metvalues_1(:,j) = s_t_1_c(:,j);
      metvalues_2(:,j) = s_t_2_c(:,j);
      gp.actstyle   = 'none';
      gp.arrowstyle = 'none';
    case 'reactions'
      actvalues_1(:,j) =  s_t_1_v(:,j);
      actvalues_2(:,j) =  s_t_2_v(:,j);
      gp.metstyle = 'none';
    otherwise, error('Unknown function option');
  end
end

jj = j;

% scale actvalues matrices to maximum values of 1
actvalues_1 = 1/max(abs(actvalues_1(:))) * actvalues_1;
actvalues_2 = 1/max(abs(actvalues_2(:))) * actvalues_2;

if length(metvalues_1), gp.metvaluesmax = max(max(metvalues_1(:)),0.001); end
if length(actvalues_1), gp.actvaluesmax = max(max(actvalues_1(:)),0.001); end

clf;
subplot('position',[0 0 1 1]);
gp.hold_on=0; gp.circle_shift = [0,0];
netgraph_concentrations(network,metvalues_1(:,1),actvalues_1(:,1),0,gp); hold on;

gp.hold_on=1; gp.circle_shift = gp.set_circle_shift;
netgraph_concentrations(network,metvalues_2(:,1),actvalues_2(:,1),0,gp);

gp_default.figure_axis    = axis;
gp_default.actprintvalues = 0;
gp_default.metprintvalues = 0;
gp_default.arrowstyle     = 'none';

gp = join_struct(gp_default,gp);

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
  gp.hold_on = 0; gp.circle_shift=[0,0];
  if length(gp.colormap_1), gp.colormap = gp.colormap_1; end
  if ~isempty(gp.background_colors)
    gp.canvas = gp.background_colors(j,:);
  end
  gp.suppress_lines = 0;
  netgraph_concentrations(network,metvalues_1(:,j),actvalues_1(:,j),text_flag,gp); hold on;
  gp.hold_on=1; gp.circle_shift=gp.set_circle_shift;
  if length(gp.colormap_2), gp.colormap = gp.colormap_2; end
  %gp.canvas = [];
  gp.suppress_lines = 1;
  netgraph_concentrations(network,metvalues_2(:,j),actvalues_2(:,j),text_flag,gp);

  if gp.timebar
    ss = network.graphics_par.squaresize;
    netgraph_movie_timebar(a, ss, (j-1)/max(jj-1,1), ss/2);
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
