% NETGRAPH_MARK_METABOLITES Mark metabolites with the mouse
%
% [choice, p] = netgraph_mark_metabolites(network, choice, draw_only)
%
% Mouse buttons: left/centre: mark/unmark a metabolite; right: finish
%
% Inputs:
%   network    network structure
%   choice     (optional) initial choice (vector of 0 and 1, one entry per metabolite)
%   draw_only  (optional, default 0) only show the current choice
%
% Outputs:
%   choice  marked metabolites (vector of 0 and 1)
%   p       graphics settings used
%
% See also: netgraph_mark_reactions, netgraph_external

function [choice,p] = netgraph_mark_metabolites(network,choice,draw_only)

if ~exist('choice','var'), choice = zeros(size(network.metabolites)); end
if length(choice)~=length(network.metabolites), error('mnt:netgraph:wrongSize', 'Wrong size of input argument'); end

p = network.graphics_par;

splitted = 0;
if isfield(p,'display_split'),   if p.display_split,     splitted =1;   end;  end

p.metvalues=0.2*ones(size(p.metnames));
p.arrowvalues=zeros(size(p.arrowvalues));
p.arrowstyle='fluxes';
p.actvalues=zeros(length(network.actions),1);
p.metvalues(find(choice))=0.5;
p.metvalues_std=zeros(length(p.metvalues),1);
p.metvaluesmax=1;
p.actvaluesmax=1;
p.metstyle='box';
p.actstyle='none';
p.actprintnames=0;

netgraph_draw(network,p); title('Left button > mark, Right button > finish')

if ~exist('draw_only','var'), draw_only = 0; end
if draw_only,  return; end

n_met = length(network.metabolites);

if splitted
  x = p.x_split(:,1:length(p.split_back_mapping));
else
  x = p.x(:,1:n_met);
end

n_met = length(network.metabolites);
cont = 1;

while cont==1
  [x_old,y_old,button] = ginput(1);
  dist = sum( (repmat([x_old;y_old],1,size(x,2))-x).^2);
  [dum,ind]=min(dist);
  index = ind(1);
  if splitted
    index = p.split_back_mapping(index);
  end
  switch button
    case 3
      cont = 0;
    otherwise
      if index <= n_met
        choice(index)=1-choice(index);
      end
  end
  p.metvalues=0.2*ones(length(p.metnames),1);
  p.metvalues(find(choice))=0.5;
  p.metvaluesmax = 1;
  p.actvaluesmax = 1;
  netgraph_draw(network,p);  title('Left button > mark, Right button > finish')
end
