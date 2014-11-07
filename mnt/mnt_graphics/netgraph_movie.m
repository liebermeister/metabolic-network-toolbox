% [M,T] = netgraph_movie(network, t ,s_t, data_type, n_frames, text_flag, options);
%
% s_t must be #variables x #timepoints !!!
%
% data_type {'concentrations','reactions','both','compute reactions'}
% default: concentrations
%
% if 'both' is chosen, then s_t must contain the matrix of concentration 
% time courses and below the matrix of reaction velocity time courses 
%
% Symbols colors for each frame can be specified in options.metabolite_colors and options.reaction_colors
%
% usage: movie(M);
%
% to save a movie as an animated gif, use movie_save(filename,M)

function [M,T] = netgraph_movie(network,t,s_t,data_type,n_frames,text_flag,options);

eval(default('options','struct','data_type','''concentrations''','n_frames','[]','text_flag','0','timebar','1'));

if ~exist('options','var'), options = struct; end

options_default = struct('flux_moves',0,'metabolite_colors',[],'reaction_colors',[]);

options = join_struct(options_default, options);

if isempty(n_frames), n_frames = 10; end

t = t-min(t);

T = max(t)*(0:1/(n_frames-1):1);

nandata = find(sum(isfinite(s_t),2)==0);
s_t(sum(isfinite(s_t),2)==0,:)=0;
s_t            = interp1(t,s_t',T,'cubic')';
s_t(nandata,:) = nan;
if size(s_t,2)==1, s_t = s_t'; end


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

[nm,nr] = size(network.N);

if isfield(network.graphics_par,'reaction_mapping'), 
  nr = max(network.graphics_par.reaction_mapping); 
  nm = size(s_t,1) - nr;
end

if isfield(network.graphics_par,'metabolite_mapping'), 
  nm = max(network.graphics_par.metabolite_mapping); 
end

metvalues = zeros(nm,length(T));
actvalues = zeros(nr,length(T));

for j = 1:length(T),
  switch data_type,    
    case 'concentrations',
      metvalues(:,j) = s_t(:,j);
      options.actstyle = 'none';
    case 'reactions',
      actvalues(:,j) =  s_t(:,j);
      options.metstyle = 'none';
    case 'both',
      metvalues(:,j) = s_t(1:nm,j);
      actvalues(:,j) = s_t(nm+1:end,j);
    case 'compute reactions',
      metvalues(:,j) = s_t(:,j);
      actvalues(:,j) = network_velocities(metvalues(:,j),network,network.kinetics);    
  end
end

jj = j;

if length(metvalues), options_default.metvaluesmax = max(max(metvalues(:)),0.001); end
if length(actvalues), options_default.actvaluesmax = max(max(actvalues(:)),0.001); end

netgraph_concentrations(network,metvalues(:,1),actvalues(:,1),0,options_default);

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
  
  netgraph_concentrations(network,metvalues(:,j),actvalues(:,j),text_flag,options);

  if timebar, 
    hold on;  
    ss = network.graphics_par.squaresize;
    plot([a(1)+ss,a(2)-ss],ss+[a(3),a(3)],'-','Color',[0.4 0.7 1]); 
    circle( a(1) + ss + (j-0.9)/(jj-0.8) *(a(2)-a(1)-2*ss),ss+a(3),ss/2,[0.2 0.5 1]);
    circle( a(1) + ss + (j-0.9)/(jj-0.8) *(a(2)-a(1)-2*ss),ss+a(3),ss/3,[0.4 0.7 1]);
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
