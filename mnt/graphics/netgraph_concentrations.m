% NETGRAPH_CONCENTRATIONS Draw a network with metabolite levels and fluxes
%
% netgraph_concentrations(network, S, J, text_flag, options)
% netgraph_concentrations(network, S, J, text_flag, 'name1', value1, ...)
% netgraph_concentrations(network, S, J, options)                    (text_flag = 0)
% netgraph_concentrations(network, S, J, 'name1', value1, ...)  (text_flag = 0)
%
% S          metabolite levels (shown by coloured circles; [] for none)
% J          fluxes (shown by arrows; [] for none)
% text_flag  show names of metabolites and reactions (default 0)
% options    graphics options (struct or name-value pairs; see 'help netgraph_draw'), plus
%              scale_arrowvalues  (default 1) scale the fluxes to a maximal absolute value of 1
%
% Example: netgraph_concentrations(network, c, v, 1, 'colorbar', 1)
%
% See also: netgraph_draw

function netgraph_concentrations(network,S,J,text_flag,varargin)

if ~exist('S','var'), S = []; end
if ~exist('J','var'), J = []; end
if ~exist('text_flag','var'), text_flag = 0; end

% options given in place of text_flag
if isstruct(text_flag) || ischar(text_flag) || isstring(text_flag)
  varargin  = [{text_flag}, varargin];
  text_flag = 0;
end

options = netgraph_options(varargin, network);

options_default = struct('actstyle','none','arrowstyle','none','arrowvalues',[],'actprintnames',0,'flag_edges',1,'arrowvaluesmax',max(abs(J)),'canvas',[],'scale_arrowvalues',1,'keep_subplot',0);

if length(J)
  options_default.arrowstyle = 'fluxes';
end

if isfield(network,'graphics_par')
  options_default = join_struct(network.graphics_par,options_default);
end

options = join_struct(options_default, options);

% (as before: text_flag = [] leaves the settings in graphics_par unchanged)
if ~isempty(text_flag) && ~text_flag
  options = join_struct(options, struct('metprintnames',0,'actprintnames',0));
end

options.actvalues = J;
options.metvalues = S;

if strcmp(options.arrowstyle,'fluxes')
  if isempty(options.arrowvalues)
    options.arrowvalues = J;
  end
end

if isfield(options,'arrowvalues')
  options.arrowvaluesmax = max(abs(options.arrowvalues));
end

if strcmp(options.actstyle, 'none')
  if ~strcmp(options.arrowstyle,'none')
    if isempty(options.arrowvalues)
      options.arrowvalues = J;
      options.arrowstyle = 'fluxes';
    end
  end
end

if options.scale_arrowvalues
  if length(options.arrowvalues) && max(abs(options.arrowvalues(:))) > 0
    options.arrowvalues = 1/max(abs(options.arrowvalues(:))) * options.arrowvalues;
  end
end

options_draw = struct('metstyle','fixed','metvalues',S,'actvalues',J,'arrowcolor',[0.85 0.1 0.1],'linecolor',[0.3 0.3 0.3]);
if ~isempty(J)
  if length(J)==1, J=J*ones(size(network.actions)); end
  options_draw = join_struct(options_draw,struct('arrowvalues',J,'arrowstyle','fluxes'));
end

options_draw = join_struct(options_draw,options);

netgraph_draw(network, options_draw);

if isempty(options_draw.canvas)
  set(gcf,'Color',[1 1 1]);
else
  set(gcf,'Color',options_draw.canvas);
end
