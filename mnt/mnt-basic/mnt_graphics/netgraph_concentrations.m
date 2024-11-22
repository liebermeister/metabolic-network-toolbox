% netgraph_concentrations(network,S,J,flag_text,options)
%
% display concentrations and (optionally) fluxes

function netgraph_concentrations(network,S,J,flag_text,options)

eval(default('S','[]','J','[]','flag_text','0','options','struct'));

opt_def = struct('actstyle','none','arrowstyle','none','arrowvalues',[],'actprintnames',0,'flag_edges',1,'arrowvaluesmax',max(abs(J)),'canvas',[],'scale_arrowvalues',1,'keep_subplot',0);%,'colormap',rb_colors);

if length(J),
  opt_def.arrowstyle = 'fluxes';
end

if isfield(network,'graphics_par'),
  opt_def = join_struct(network.graphics_par,opt_def);
end

options = join_struct(opt_def, options);

if ~flag_text,
  options = join_struct(options, struct('metprintnames',0,'actprintnames',0));
end

options.actvalues = J;
options.metvalues = S;

if strcmp(options.arrowstyle,'fluxes')
  if strcmp(options.arrowvalues,[]),
    options.arrowvalues = J;
  end
end

if isfield(options,'arrowvalues'),
  options.arrowvaluesmax = max(abs(options.arrowvalues));
  %options.arrowstyle = 'fluxes';
  %options.actstyle   = 'fluxes';
end

if length(J), 
  %options.actstyle = 'fixed';
end

if strcmp(options.actstyle, 'none'),
  if ~strcmp(options.arrowstyle,'none'),
    if isempty(options.arrowvalues),
      options.arrowvalues = J;
      options.arrowstyle = 'fluxes';
    end
  end
end

if options.scale_arrowvalues,
  if length(options.arrowvalues),
    options.arrowvalues = 1/nanmax(abs(options.arrowvalues(:))) * options.arrowvalues;
  end
end

opt = struct('metstyle','fixed','metvalues',S,'actvalues',J,'arrowcolor',[1 0 0],'linecolor',[0 0 0]);
if isempty(J),
%  opt = join_struct(opt,struct('arrowstyle','none'));
else
  if length(J)==1, J=J*ones(size(network.actions)); end
  opt = join_struct(opt,struct('arrowvalues',J,'arrowstyle','fluxes'));
end

opt = join_struct(opt,options);
if opt.keep_subplot==0,
  subplot('Position', [0 0 1 0.90]);
end

netgraph_draw(network, opt);

if isempty(opt.canvas),
  set(gcf,'Color',[1 1 1]);
else
  set(gcf,'Color',opt.canvas);
end
