% NETGRAPH_CONCENTRATIONS Draw a network with metabolite levels and fluxes
%
% netgraph_concentrations(network, S, J, text_flag, gp)
% netgraph_concentrations(network, S, J, text_flag, 'name1', value1, ...)
% netgraph_concentrations(network, S, J, gp)                    (text_flag = 0)
% netgraph_concentrations(network, S, J, 'name1', value1, ...)  (text_flag = 0)
%
% S          metabolite levels (shown by coloured circles; [] for none)
% J          fluxes (shown by arrows; [] for none)
% text_flag  show names of metabolites and reactions (default 0)
% gp         graphics options (struct or name-value pairs; see 'help netgraph_draw'), plus
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

gp = netgraph_options(varargin, network);

gp_default = struct('actstyle','none','arrowstyle','none','arrowvalues',[],'actprintnames',0,'flag_edges',1,'arrowvaluesmax',max(abs(J)),'canvas',[],'scale_arrowvalues',1,'keep_subplot',0);

if length(J)
  gp_default.arrowstyle = 'fluxes';
end

if isfield(network,'graphics_par')
  gp_default = join_struct(network.graphics_par,gp_default);
end

gp = join_struct(gp_default, gp);

% (as before: text_flag = [] leaves the settings in graphics_par unchanged)
if ~isempty(text_flag) && ~text_flag
  gp = join_struct(gp, struct('metprintnames',0,'actprintnames',0));
end

gp.actvalues = J;
gp.metvalues = S;

if strcmp(gp.arrowstyle,'fluxes')
  if isempty(gp.arrowvalues)
    gp.arrowvalues = J;
  end
end

if isfield(gp,'arrowvalues')
  gp.arrowvaluesmax = max(abs(gp.arrowvalues));
end

if strcmp(gp.actstyle, 'none')
  if ~strcmp(gp.arrowstyle,'none')
    if isempty(gp.arrowvalues)
      gp.arrowvalues = J;
      gp.arrowstyle = 'fluxes';
    end
  end
end

if gp.scale_arrowvalues
  if length(gp.arrowvalues) && max(abs(gp.arrowvalues(:))) > 0
    gp.arrowvalues = 1/max(abs(gp.arrowvalues(:))) * gp.arrowvalues;
  end
end

gp_draw = struct('metstyle','fixed','metvalues',S,'actvalues',J,'arrowcolor',[0.85 0.1 0.1],'linecolor',[0.3 0.3 0.3]);
if ~isempty(J)
  if length(J)==1, J=J*ones(size(network.actions)); end
  gp_draw = join_struct(gp_draw,struct('arrowvalues',J,'arrowstyle','fluxes'));
end

gp_draw = join_struct(gp_draw,gp);

netgraph_draw(network, gp_draw);

if isempty(gp_draw.canvas)
  set(gcf,'Color',[1 1 1]);
else
  set(gcf,'Color',gp_draw.canvas);
end
