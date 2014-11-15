
% ----------------------------------------------------------------------------
% graphics

if isempty(me_options.arrowvaluesmax),
  arrowvaluesmax = min(max(abs(res.v_proj)), 5* median(abs(res.v_proj)));
else
  arrowvaluesmax = me_options.arrowvaluesmax;
end

network_aug_CoHid.graphics_par.omitreactions   = me_options.omitreactions;
network_aug_CoSplit.graphics_par.omitreactions = me_options.omitreactions;


nk = length(kinetic_models);
figure(1); netgraph_concentrations(network_aug_CoHid,res.covered_metabolites + 0.5*res.shared_metabolites, res.covered_reactions + 0.5*res.shared_reactions,1,struct('actprintnames',0, 'actstyle','fixed','arrowstyle','none','text_offset',me_options.text_offset,'showsign',0,'metvaluesmax',nk,'actvaluesmax',nk,'metvaluesmin',0,'actvaluesmin',0,'colormap',[1 1 1; copper],'single_arrow',1));

for it = 1:length(kinetic_models),
  my_external = nan * ones(length(network.metabolites),1);
  my_external(res.mapping_metabolites{it}) = kinetic_models{it}.external;
  my_prod = nan * ones(length(network.metabolites),1);
  my_prod(res.mapping_metabolites{it}) = kinetic_models{it}.N * res.v_stat{it};
  my_v = res.collect_v(:,it);
  my_c = res.collect_c(:,it);  
  figure(100+it); netgraph_concentrations(network_aug_CoHid,  my_external,my_v,1,struct('arrowsize',me_options.arrowsize,'actstyle','none','text_offset',me_options.text_offset,'single_arrow',1));
  title(sprintf('Embedded model %d: External metabolites and fluxes',it))
  figure(110+it); netgraph_concentrations(network_aug_CoHid,  my_prod,my_v,1,struct('arrowsize',me_options.arrowsize,'actstyle','none','text_offset',me_options.text_offset,'single_arrow',1));
  title(sprintf('Embedded model %d: Net production and fluxes',it))
  figure(120+it); netgraph_concentrations(network_aug_CoHid,  double(isfinite(my_c)),double(isfinite(my_v)),1,struct('arrowstyle','none','text_offset',me_options.text_offset,'single_arrow',1));
  title(sprintf('Embedded model %d: Subnetwork covered',it))
end

figure(2); netgraph_concentrations(network_aug_CoHid,network_aug.external,res.v_mean,1,struct('arrowvaluesmax',arrowvaluesmax,'arrowsize',me_options.arrowsize,'squaresize',me_options.squaresize,'actstyle','none','text_offset',me_options.text_offset,'single_arrow',1));

figure(3); netgraph_concentrations(network_aug_CoHid,network_aug.external,res.v_proj,1,struct('arrowvaluesmax',arrowvaluesmax,'arrowsize',me_options.arrowsize,'squaresize',me_options.squaresize,'actstyle','none','text_offset',me_options.text_offset,'single_arrow',1));

figure(203); set(gcf, 'Position', [120 150 700 850]);
M = netgraph_flux_movie(network_aug_CoHid,network_aug.external,res.v_proj,1,struct('arrowvaluesmax',arrowvaluesmax,'arrowsize',me_options.arrowsize,'squaresize',me_options.squaresize,'actstyle','none','text_offset',me_options.text_offset,'metprintnames',0));
%%movie(M,10)

figure(204); set(gcf, 'Position', [120 150 700 850]);
[M2,T] = netgraph_movie(network_aug_CoHid,simulation.t,log(diag(1./simulation.C(:,end))*simulation.C),'concentrations',40,0);
 

figure(13); netgraph_concentrations(network_aug_CoHid,res.c_fix,[],1,struct('text_offset',me_options.text_offset,'single_arrow',1));

figure(14); netgraph_concentrations(network_aug_CoHid,[],res.dmu_fix,1,struct('text_offset',me_options.text_offset,'arrow_stoichiometries',0,'single_arrow',1));

for it = 1:length(kinetic_models),
  figure(200+it); netgraph_concentrations(network_aug_CoHid,res.collect_mu(:,it),[],1,struct('text_offset',me_options.text_offset,'arrow_stoichiometries',0,'single_arrow',1));
end

figure(4); netgraph_concentrations(network_aug_CoHid,res.mu-min(res.mu),res.A,1,struct('arrowvaluesmax',max(abs(res.A)),'arrowsize',me_options.arrowsize,'squaresize',me_options.squaresize,'actstyle','none','text_offset',me_options.text_offset,'colormap',flipud(rb_colors),'single_arrow',1));

figure(5); netgraph_concentrations(network_aug_CoHid,res.c_combined,res.v_combined,1,struct('arrowvaluesmax',arrowvaluesmax,'arrowsize',me_options.arrowsize,'squaresize',me_options.squaresize,'actstyle','none','text_offset',me_options.text_offset,'single_arrow',1));

figure(6); 
h = plot(simulation.t,log10(simulation.C)); 
line_colors(h,'rby_colors');
legend(network.metabolites,'FontSize',me_options.legendFontsize);
xlabel('Time'); ylabel('Log10 concentration');

figure(7);
if length(simulation_flat.t),
  h = plot(simulation_flat.t,log10(simulation_flat.C)); 
  line_colors(h,'rby_colors');
  legend(network.metabolites,'FontSize',me_options.legendFontsize);
  xlabel('Time'); ylabel('Log10 concentration');
end


figure(8);
h = plot(simulation_local.t,log10(simulation_local.C)); 
line_colors(h,'rby_colors');
legend(kinetic_models{1}.metabolites,'FontSize',me_options.legendFontsize);
xlabel('Time'); ylabel('Log10 concentration');

figure(9);
h = plot(simulation.t,log10(simulation.C(res.mapping_metabolites{1},:))); 
line_colors(h,'rby_colors');
legend(network_aug.metabolites(res.mapping_metabolites{1}),'FontSize',me_options.legendFontsize);
xlabel('Time'); ylabel('Log10 concentration');


% ----------------------------------------------------------------------------
% print graphics

if flag_save_graphics,

cd(filenames.psfile_dir);

print embedding_kinetic_mapping.eps      -f1 -depsc
print embedding_kinetic_flux_kinetic.eps -f2 -depsc
print embedding_kinetic_flux_network.eps -f3 -depsc
print embedding_kinetic_thermo.eps       -f4 -depsc
print embedding_kinetic_result_state.eps -f5 -depsc
print embedding_kinetic_integration.eps  -f6 -depsc
print embedding_kinetic_integration_flat.eps  -f7 -depsc
print embedding_kinetic_integration_isolated.eps  -f8 -depsc
print embedding_kinetic_integration_embedded.eps  -f9 -depsc

for it = 1:length(kinetic_models),
  print(['embedding_kinetic_model_' num2str(it) '_fluxes.eps'],['-f' num2str(100+it)], '-depsc');
  print(['embedding_kinetic_model_' num2str(it) '_mapping.eps'],['-f' num2str(120+it)], '-depsc');
end

display(sprintf('Saving graphics to directory %s',filenames.psfile_dir));

cd(filenames.psfile_dir);
movie_file = 'embedding_kinetic_movie';
movie_save(movie_file,M)
display(sprintf('Saving flux movie to file %s/%s.gif',filenames.psfile_dir,movie_file));

movie_file = 'embedding_kinetic_movie_perturbation';
movie_save(movie_file,M2)
display(sprintf('Saving flux movie to file %s/%s.gif',filenames.psfile_dir,movie_file));

end

figure(1); title('Kinetic models embedded');
figure(2); title('Flux projection: predefined fluxes and external metabolites');
figure(3); title('Flux projection: stationary fluxes and external metabolites');
figure(13); title('Concentrations taken from kinetic models');
figure(14); title('Fixed delta mu');
for it = 1:length(kinetic_models),
  figure(200+it);  title('Fixed chemical potentials mu');
end
figure(4); title('Chemical potentials (shifted to min=0) and reaction affinities');
figure(5); title('Combined network: steady state concentrations and fluxes');
figure(6); title('Combined network: simulation');
figure(7); title('Combined network: simulation of flattened model');
figure(8); title('Simulation of isolated kinetic model');
figure(9); title('Simulation of embedded kinetic model');
