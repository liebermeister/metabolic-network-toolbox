% [M,T] = netgraph_movie(network, t, s_t_1, s_t_2, data_type, n_frames, text_flag, options);
%
% Movie displaying data on network; two data sets are shown at teh same time 
% (for instance, to compare experimental and simulated time series)

function [M,T] = netgraph_movie(network,t,s_t_1,s_t_2,data_type,n_frames,text_flag,options);

eval(default('options','struct','data_type','''concentrations''','n_frames','[]','text_flag','0','timebar','1'));

if ~exist('options','var'), options = struct; end

options_default = struct('flux_moves',0,'metabolite_colors',[],'reaction_colors',[]);

options = join_struct(options_default, options);

if isempty(n_frames), n_frames = 10; end

t = t-min(t);

T = max(t)*(0:1/(n_frames-1):1);


nandata = find(sum(isfinite(s_t_1),2)==0);
s_t_1(sum(isfinite(s_t_1),2)==0,:)=0;
s_t_1            = interp1(t,s_t_1',T,'cubic')';
s_t_1(nandata,:) = nan;
if size(s_t_1,2)==1, s_t_1 = s_t_1'; end


nandata = find(sum(isfinite(s_t_2),2)==0);
s_t_2(sum(isfinite(s_t_2),2)==0,:)=0;
s_t_2            = interp1(t,s_t_2',T,'cubic')';
s_t_2(nandata,:) = nan;
if size(s_t_2,2)==1, s_t_2 = s_t_2'; end


if length(options.metabolite_colors), 
  options.metabolite_colors = interp1(t,options.metabolite_colors,T,'cubic');
  options.metabolite_colors(options.metabolite_colors<0) = 0;
  options.metabolite_colors(options.metabolite_colors>1) = 1;
end

if length(options.reaction_colors), 
  options.reaction_colors = interp1(t,options.reaction_colors,T,'cubic');
  options.reaction_colors(options.reaction_colors<0) = 0;
  options.reaction_colors(options.reaction_colors>1) = 1;
end

switch data_type,    
  case 'concentrations',
    s_t_1_c = s_t_1;
    s_t_1_v = [];
    s_t_2_c = s_t_2;
    s_t_2_v = [];
    nm    = size(s_t_1,1);
    nr    = 0;
  case 'reactions',
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

for j = 1:length(T),
  switch data_type,    
    case 'concentrations',
      metvalues_1(:,j) = s_t_1_c(:,j);
      metvalues_2(:,j) = s_t_2_c(:,j);
      options.actstyle = 'none';
    case 'reactions',
      actvalues_1(:,j) =  s_t_1_v(:,j);
      actvalues_2(:,j) =  s_t_2_v(:,j);
      options.metstyle = 'none';
   otherwise, error('Unknown function option'); 
  end
end

jj = j;

if length(metvalues_1), options_default.metvaluesmax = max(max(metvalues_1(:)),0.001); end
if length(actvalues_2), options_default.actvaluesmax = max(max(actvalues_1(:)),0.001); end

clf; 
subplot('position',[0 0 1 1]);
options.hold_on=0; options.circle_shift=[0,0];
netgraph_concentrations(network,metvalues_1(:,1),actvalues_1(:,1),0,options_default); hold on;
options.hold_on=1; options.circle_shift=[0.01,0.01];
netgraph_concentrations(network,metvalues_2(:,1),actvalues_2(:,1),0,options_default);

options_default.figure_axis    = axis;
options_default.actprintvalues = 0;
options_default.metprintvalues = 0; 
options_default.arrowstyle     = 'none';

options = join_struct(options_default,options);

if isfield(options,'timebar'), timebar = options.timebar; end
a = axis;
clf;

for j=1:jj,

  set(gcf,'Color',[1 1 1]);

  if options.flux_moves,
    arrow_shift         = mod(options.flux_moves*j,jj)/jj;
    options.arrow_shift = arrow_shift;
  end

  if length(options.metabolite_colors),
    options.metcolors = squeeze(options.metabolite_colors(j,:,:));
    if size( options.metcolors,2) == 1, options.metcolors = options.metcolors'; end
  end

  if length(options.reaction_colors),
    options.actcolors = squeeze(options.reaction_colors(j,:,:));
    if size( options.actcolors,2) == 1, options.actcolors = options.actcolors'; end
  end
  
  clf; 
  subplot('position',[0 0 1 1]);
  options.hold_on=0; options.circle_shift=[0,0];
  netgraph_concentrations(network,metvalues_1(:,j),actvalues_1(:,j),text_flag,options); hold on;
  options.hold_on=1; options.circle_shift=[0.01,0.01];
  netgraph_concentrations(network,metvalues_2(:,j),actvalues_2(:,j),text_flag,options);
  if timebar, 
    hold on;  
    ss = network.graphics_par.squaresize;
    plot([a(1)+ss,a(2)-ss],ss+[a(3),a(3)],'-','Color',[0 0 0]); 
    circle( a(1) + ss + (j-0.9)/(jj-0.8) *(a(2)-a(1)-2*ss),ss+a(3),ss/2,[0 0 0]);
    circle( a(1) + ss + (j-0.9)/(jj-0.8) *(a(2)-a(1)-2*ss),ss+a(3),ss/6,[1 1 1]);
  end

  axis(a);

  M(j) = getframe(gca);

  if isfield(options,'prefix'),
    if j ==1,             print([options.prefix '_start.eps'],'-f1','-depsc'); end
    if j == ceil(0.8*jj), print([options.prefix '_end.eps'],'-f1','-depsc');   end
  end

end

if isfield(options,'rest_frames'),
  M = [ repmat(M(1),1,options.rest_frames), M, repmat(M(end),1,options.rest_frames)];
end
