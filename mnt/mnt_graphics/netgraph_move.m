%network = netgraph_move(network,flag_draw_details,options,picture,filename_in,filename_out)
%
%move nodes of metabolic network graph with the mouse
%
% left button:     first click: choose node
%                 second click: choose new position
% center button:   relax positions
% right button:    quit
%
% flag_draw_details (optional): 
%  'none': only reaction lines
%  'regulation': also regulation lines
%  'text': also node names
%  'all':  also regulation lines + node names
%
% options: subfield 'initial_relax' leads to initial relaxation_steps
% options: subfield 'manual' = 0 -> only automatic relaxation 

function network = netgraph_move(network,flag_draw_details,options,picture,filename_in,filename_out)

if ~exist('flag_draw_details','var'), flag_draw_details = 'text';  end;
if ~exist('picture','var'), picture = []; end 
if ~exist('options','var'),  options = struct; end
if ~isfield(options,'manual'), options.manual = 1; end

if isempty(flag_draw_details), flag_draw_details = 'none';  end;

if exist('filename_in','var'),
  file = fopen(filename_in);
  it = 1;
  clear x names
  while ~feof(file),
    lin   = fgetl(file);
    start = regexp(lin,'[^\s]+');
    stop  = regexp(lin,'[\s]+');
    names{it} = lin(start(1):stop(1)-1);
    x(1,it) = eval(lin(start(2):stop(2)-1)); 
    x(2,it) = eval(lin(start(3):stop(3)-1));
    fixed(it) = eval(lin(start(4):end));
    it = it+1;
  end
  fclose(file);
  names= names';
  ll = label_names(names,network.metabolites);
  network.graphics_par.x(:,ll) = x;
  network.graphics_par.fixed(ll) = fixed;
end

switch flag_draw_details,
  case 'none', opt = struct('metstyle','none','actstyle','none','metprintnames',0,'actprintnames',0,'metprintvalues',0,'actprintvalues',0);
  case 'regulation', opt = struct('metstyle','none','actstyle','none','metprintnames',0,'actprintnames',0,'metprintvalues',0,'actprintvalues',0,'show_regulation',1); 
  case {'text','all'}, opt = struct('metstyle','none','actstyle','none'); 
end

opt = join_struct(opt,options);

if ~isfield(network,'graphics_par'),  network=netgraph_make_graph(network); end

% ------------------------------------------------------------
% layout parameters

if ~isfield(opt,'frame'),
  p = network.graphics_par;
  if sum(sum(p.x ~= 0 )) == 0, opt.frame = [0 1 0 1]; end 
  opt.frame = [1.1 * min(p.x(1,:))-0.1* max(p.x(1,:)),...
      -0.1 * min(p.x(1,:))+1.1* max(p.x(1,:)),...
      1.1 * min(p.x(2,:))-0.1* max(p.x(2,:)),...
      -0.1 * min(p.x(2,:))+1.1* max(p.x(2,:))];
% frame of the model
end

a = opt.frame;

nm        = length(network.graphics_par.metnames);
size_unit = sqrt((a(2)-a(1))*(a(4)-a(3))) * sqrt(1/nm);  % typical
                                                         % distance
                                                         % between equally spaced dots
force_unit = 0.1;

gridsize = 0.01  * size_unit;   % distance between discrete x values
k1       = 0.005 * force_unit;  % overall repulsion strength
k2       = 0.1   * size_unit;   % overall repulsion scale
k3       = 0.5   * force_unit;  % neighbour attraction strength
k4       = 0.5   * size_unit;   % neighbour attraction scale

p   = network.graphics_par;

splitted = 0; if isfield(p,'display_split'), splitted = p.display_split;   end; 

if splitted,
  if ~isfield(network.graphics_par,'fixed_split'), 
    p.fixed_split = zeros(size(network.graphics_par.x_split,2),1); 
  end
  fixed = p.fixed_split;
  N   = p.N_split;
  x   = p.x_split;
  m   = p.m_split;
  n_met = length(p.split_back_mapping);
else,
  if ~isfield(p,'fixed'), p.fixed = zeros(length(network.metabolites)+length(network.actions),1); end
  fixed = p.fixed;
  N  = p.N;
  x   = p.x;
  m   = p.m;
  n_met = length(p.metnames);
end


switch flag_draw_details,
  case {'regulation','all'},
    if ~isfield(p,'regulation_matrix'); flag_draw_details = 'none'; end
    p.show_regulation = 1; network.graphics_par.show_regulation = 1; 
end

p.arrowstyle='none';

ind_hub_metabolites = find(sum(double(network.N~=0),2)>5);

% -----------------------------------------------------------------
% first relaxation and display

names=[p.metnames; p.actnames];
loose = find(1-fixed);

if isfield(opt,'initial_relax'),
  for it = 1:opt.initial_relax, x = netgraph_relaxation_step(x,loose,m,k1,k2,k3,k4,gridsize,ind_hub_metabolites); end
end

display_graphics(picture,network,p,opt,x,fixed);

% -----------------------------------------------------------------
% manual relaxation

if opt.manual ~=0, cont = 1; else, cont = 0; end

[Ni,Nj]=ind2sub(size(N),find((N~=0)));
Nj = Nj+n_met;

%edgecenters = 0.5*(x(:,Ni)+x(:,Nj));
%hold on;      %     plot(edgecenters(1,:),edgecenters(2,:),'g.');
%fixed_ind = find(fixed);
%plot(x(1,fixed_ind),x(2,fixed_ind),'k.');
%hold off;
%loose = find(1-fixed);

% m is the adjacency matrix

while cont==1,

  button = 0;  
  [x_old,y_old,button] = ginput(1);

  while isempty(button),
    [x_old,y_old,button] = ginput(1);
  end

  switch button
    
    case 1,
      dist    = sum( (repmat([x_old;y_old],1,size(x,2))-x).^2);
      [dum,i] = min(dist);
      index   = i(1);
      title('Left > place node. Center > place and fix node');
      [x_new,y_new,button] = ginput(1);
            x(:,index) = [x_new;y_new];
%      factors = exp(-db(index,:)/p.lambda);
%      x(:,loose) = x(:,loose) + [x_new-x_old; y_new-y_old] * factors(:,loose);
      if  button-1 ~= fixed(index),
	fixed(index) = button-1;
	loose = find(1-fixed);
      end
      
    case 2,  for it = 1:20, x = netgraph_relaxation_step(x,loose,m,k1,k2,k3,k4,gridsize,ind_hub_metabolites); end

    case 3,   cont = 0; 
  
  end

%  edgecenters = 0.5*(x(:,Ni)+x(:,Nj));
  
  if splitted,     p.x_split = x;
  else,            p.x = x;
  end
  
  network.graphics_par=p;

  display_graphics(picture,network,p,opt,x,fixed);  
  
end

if exist('filename_out','var'), netgraph_save_positions(network,filename_out); end
if opt.manual, title('New positions confirmed'); end

p.x      = x; 
p.fixed = fixed; 

if splitted,
  p.x_split     = x; 
  p.fixed_split = fixed; 
end

network.graphics_par = p;

% --------------------------------------------------------------------------------

function display_graphics(picture,network,p,opt,x,fixed)

if opt.manual,

network.graphics_par = p;
network.graphics_par.x = x;

if exist('picture','var'), 
  image(picture); 
  hold on;
end;

netgraph_draw(network,opt); hold on;
title('Left > move node. Center > relax network. Right > finish');
drawnow; 

fixed_ind = find(fixed);
plot(x(1,fixed_ind),x(2,fixed_ind),'k.'); 
hold off;
if isfield(opt,'Axis'), 
  axis(opt.Axis); 
else,
  axis(opt.frame);
end
axis tight
end

% --------------------------------------------------------------------------------

function  x = netgraph_relaxation_step(x,loose,m,k1,k2,k3,k4,gridsize,ind_hub_metabolites)

x(:,loose) = x(:,loose) + 0.001*randn(2,length(loose));

% spring forces:
f2 = (( m(loose,:)*x')- repmat(sum(m(loose,:),2),1,2).*x(:,loose)')';

% nonlinear springs

f2_norm = sqrt(sum(f2.^2))/k4;
f2      = f2 .*repmat( 0.5*f2_norm+1./(0.2+f2_norm),2,1);

% repelling forces
thresh = 10*k2;
d1 = abs( repmat(x(1,:),length(loose),1)-repmat(x(1,loose)',1,size(x,2))) < thresh;
d2 = abs( repmat(x(2,:),length(loose),1)-repmat(x(2,loose)',1,size(x,2))) < thresh;
%d_non_neighb = 
d3 = d1.*d2.*(1-m(loose,:));
d3(:,1:size(d3,1)) = d3(:,1:size(d3,1)) - speye( size(d3,1));
has_neighbours = find(sum(d3,2));

f0 = zeros(2,size(d3,1));      
%      for it=1:1,
%	edgecenters = 0.5*(x(:,Ni)+x(:,Nj));

for ittt=1:length(has_neighbours);
  itt = has_neighbours(ittt);
  neighbours = find(d3(itt,:));
  distances = repmat(x(:,loose(itt)),1,length(neighbours)) - x(:,neighbours);
  f1 = sum(distances.^2)+10^-10;
  f1 = exp(-f1/(2*k2^2))./sqrt(f1);
  f0(:,itt) = sum(distances.*repmat(f1,2,1),2);	  
end
%      end

%f2(:,ind_hub_metabolites) = 0;

if length(loose),
  force_loose = 0*k3 * f2 + k1 * f0;
  x(:,loose) = x(:,loose) + force_loose;
  x(:,loose) = x(:,loose)-repmat(min(x(:,loose)')',1,size(x(:,loose),2));
%  x(:,loose) = diag(1./(max(x(:,loose)')'))*x(:,loose);
%  x(:,loose) = gridsize * round(x(:,loose)/gridsize);
end

 x(1,:) = x(1,:)-min(x(1,:)); x(1,:) = x(1,:)/max(x(1,:));
 x(2,:) = x(2,:)-min(x(2,:)); x(2,:) = x(2,:)/max(x(2,:));
