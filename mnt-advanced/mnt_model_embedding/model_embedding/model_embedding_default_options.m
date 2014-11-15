function me_options = model_embedding_default_options()

me_options.id                   = [];
me_options.fba_constraints      = [];
me_options.make_kinetic_models_stationary = 1;
me_options.cofactors            = [];  % for graphics 
me_options.position_file        = [];  % for graphics 
me_options.v_network            = [];
me_options.arrowsize            = 0.02;  % for graphics 
me_options.squaresize           = 0.02;  % for graphics 
me_options.text_offset          = [0.02,0.01];  % for graphics 
me_options.arrowvaluesmax       = [];  % for graphics 
me_options.legendFontsize       = 8;  % for graphics 
me_options.fontsize             = 8;  % for graphics 
me_options.verbose              = 0;
me_options.omitreactions        = {};  % for graphics 
