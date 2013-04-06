% [choice,p] = netgraph_mark_reactions(network,choice,draw_only)

function [choice,p] = netgraph_mark_reactions(network,choice,draw_only)

% left/centre button: add/remove reversible
% right button quit

x = network.graphics_par.x;
 
  p.metstyle='none';
  p.arrowstyle='none';
  p.actstyle='box';

  p.metvalues=zeros(length(network.metabolites),1);
  p.metvalues_std=zeros(length(network.metabolites),1);
  p.actvalues=0.2*ones(length(network.actions),1);
  p.actvalues(find(choice))=0.5;
  p.actvalues_std=zeros(length(network.actions),1);
  p.metvaluesmax=1;
  p.actvaluesmax=1;
  p.arrowstyle='directions';
  netgraph_draw(network,p);
  
if ~exist('draw_only','var'), draw_only = 0; end
if draw_only,  return; end
  title('Left button > mark, Right button > finish')


  n_met=length(network.metabolites);
  cont = 1;
  
  while cont==1,
    [x_old,y_old,button] = ginput(1);
    dist= sum( (repmat([x_old;y_old],1,size(x,2))-x).^2);
    [dum,i]=min(dist);
    index=i(1);
    switch button
      case 3,
     cont = 0; 
      otherwise
	if index>n_met, choice(index-n_met)=1-choice(index-n_met); end
    end
 p.actvalues=0.2*ones(length(network.actions),1);
 p.actvalues(find(choice))=0.5;
 netgraph_draw(network,p);
 title('Left button > mark, Right button > finish')
end
