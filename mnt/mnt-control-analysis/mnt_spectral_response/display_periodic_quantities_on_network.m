function display_periodic_quantities_on_network(network_CoSplit, network, c_osc, v_osc, c_ss, v_ss, squaresize, circle_size, fontsize, fignum, metabolite_color, reaction_color);

eval(default('metabolite_color','[]','reaction_color','[0.6 0.6 1]'));
if isempty(metabolite_color), metabolite_color = [0 0 1]; end
if isempty(reaction_color), reaction_color = [0.7 0.4 0.1]; end

% phase=0 -> vertical pointer

clock_style = 1;

[nm,nr]  = size(network_CoSplit.N);

if length(fignum), figure(fignum(1)); clf; end

if length(c_osc)
  mm = network_CoSplit.graphics_par.metabolite_mapping;
  mC = circle_size / max(abs(c_osc(mm)));
  netgraph_concentrations(network_CoSplit,[],sign(v_ss),1,struct('actstyle','none','arrowstyle','fluxes','arrowvalues',sign(v_ss),'single_arrow',1,'squaresize',squaresize,'metcolstyle','fixed','metcol',metabolite_color,'FontSize',fontsize,'hold_on',1,'keep_subplot',1,'metprintnames',0,'actprintnames',0)); hold on;
  my_feather(mC*c_osc(mm),0.1,network_CoSplit.graphics_par.x(:,1:nm),[],metabolite_color,struct('show_circle',1,'clock_style',clock_style));  hold on;
  hold off
  if length(fignum), figure(fignum(2)); clf; end
end

if length(v_osc)
  rm = network_CoSplit.graphics_par.reaction_mapping;
  mV = circle_size / max(abs(v_osc(rm)));
  netgraph_concentrations(network_CoSplit, ones(size(c_ss)),[],1,struct('actstyle','none','metstyle','box','arrowstyle','none','arrowvalues',v_ss,'single_arrow',1,'squaresize',squaresize,'metcolstyle','fixed','metcol',metabolite_color,'FontSize',fontsize,'hold_on',1,'keep_subplot',1,'metprintnames',0,'actprintnames',0)); hold on;

  my_feather(mV * v_osc(rm),0.1,network_CoSplit.graphics_par.x(:,nm+[1:nr]),[],reaction_color,struct('show_circle',1,'clock_style',clock_style));

  hold off
  axis equal; axis tight; noticks
end
