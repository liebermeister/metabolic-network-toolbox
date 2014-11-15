function me_options_default = me_default_options(nk)

% nk: number of kinetic models to be embedded
%
% me_options.c_init               
% me_options.Tmax                 
% me_options.set_external         
% me_options.set_internal         
% me_options.reaction_fixed_signs 
% me_options.fixed_signs          
% me_options.id                   
% me_options.fba_constraints      
% me_options.verbose              
% me_options.position_file        
% me_options.cofactors            
% me_options.v_network            
% me_options.arrowsize            
% me_options.squaresize           
% me_options.text_offset          
% me_options.arrowvaluesmax
% me_options.legendFontsize
% me_options.make_kinetic_models_stationary


% --------------------------------------------------------------------------
% default values for me_options

me_options_default.c_init               = [];
me_options_default.Tmax                 = 5; 
me_options_default.set_external         = {};
me_options_default.set_internal         = {};
me_options_default.reaction_fixed_signs = {};
me_options_default.fixed_signs          = [];
me_options_default.id                   = [];
me_options_default.fba_constraints      = [];
me_options_default.verbose              = 0;
me_options_default.position_file        = [];
me_options_default.cofactors            = {};
me_options_default.v_network            = [];
me_options_default.arrowsize            = 0.02;
me_options_default.squaresize           = 0.02;
me_options_default.text_offset          = [0.01,0.005];
me_options_default.arrowvaluesmax       = []; 
% maximal flux (shown in graphics)
me_options_default.legendFontsize       = 4;
me_options_default.make_kinetic_models_stationary = 0;

id_default.metabolites_network        = 'metabolites';
id_default.reactions_network          = 'actions';
id_default.metabolites_kinetic_models = repmat({'metabolites'},nk,1);
id_default.reactions_kinetic_models   = repmat({'actions'},nk,1);

me_options_default.id = id_default;
