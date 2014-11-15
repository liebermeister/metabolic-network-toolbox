% --------------------------------------------------------------------------
% Run model embedding
% Small example 1 (one kinetic model, kinetics type 'kinetic_strings')

model_name = 'embedding_example_1_sbml';

filenames = filenames_embedding(model_name);

load(filenames.network_file);

eval(default('id','struct','fba_constraints','struct'));

c_init               = [0 0 0 0 1]';
Tmax                 = 5; 
me_options.id                   = id;
me_options.make_kinetic_models_stationary = 0;

me_options.arrowsize            = 0.1;
me_options.squaresize           = 0.1;
me_options.text_offset          = [0.05,0.02];
me_options.arrowvaluesmax       = 1.2;
me_options.legendFontsize       = 8;

% ------------------------------------------------------------
% Run model embedding

me_options.verbose = 0;
me_options.position_file = filenames.table_positions;

[network,network_CoHid,network_combined, res] = model_embedding(kinetic_models, network, network_CoHid, me_options );

% run simulation
% c_init = c_combined; c_init(1) = 5 * c_init(1);

[simulation.t, simulation.C] = network_integrate(network_combined, c_init, Tmax);

network_flat = model_embedding_kinetics_combined2flat(network_combined);

[simulation_flat.t, simulation_flat.C] = network_integrate(network_flat, c_init, Tmax);

flag_save_graphics = 0;

model_embedding_graphics

%network_sbml_export(network_flat,0,'Combined model',[filenames.network_dir '/' model_name '_flat.xml']);
  
