% netgraph_concentrations(network,S,J,flag_text,options)
%
% display concentrations and (optionally) fluxes

function netgraph_concentrations(network,S,J,flag_text,options)

if ~exist('S','var'),S=[]; end

if ~exist('J','var'), J = []; end

if ~exist('flag_text','var'), flag_text=0; end

opt_def = struct('actprintnames',0,'flag_edges',1);%,'colormap',rb_colors);
eval(default('options','struct'));
options = join_struct(opt_def,options);

if flag_text,
  if ~isempty(J),
    if length(J)==1,J=J*ones(size(network.actions)); end
    opt = struct('metstyle','fixed','metvalues',S,'actvalues',J,'arrowvalues',J,'arrowstyle','fluxes','arrowcolor',[0.7 0.7 0.7]);
  else
    opt = struct('metstyle','fixed','metvalues',S,'actvalues',J,'arrowstyle','none','arrowcolor',[0.7 0.7 0.7]);
  end
else,
  if ~isempty(J),
    if length(J)==1,J=J*ones(size(network.actions)); end
    opt = struct('metstyle','fixed','metvalues',S,'actvalues',J,'arrowvalues',J,'arrowstyle','fluxes','metprintnames',0,'actprintnames',0,'arrowcolor',[0.7 0.7 0.7]);
  else
    opt = struct('metstyle','fixed','metvalues',S,'arrowstyle','none','metprintnames',0,'actprintnames',0,'arrowcolor',[0.7 0.7 0.7]);
  end
end

opt = join_struct(opt,options);

netgraph_draw(network,opt);
set(gcf,'Color',[1 1 1]);