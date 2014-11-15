function [network_combined, res,network_aug,network_aug_CoHid,network_aug_CoSplit] = model_embedding(kinetic_models, network, network_CoHid, me_options)

%[network,network_CoHid,network_combined, res] = model_embedding(kinetic_models, network, network_CoHid, me_options)
%
% Inputs
%   kinetic_models
%   network
%   network_CoHid
%   me_options
%     me_options.v_network
%     me_options.id   structure describing mappings between the models. 
%                                     Default:
%       me_options.id.metabolites_network        = 'metabolites';
%       me_options.id.reactions_network          = 'actions';
%       me_options.id.metabolites_kinetic_models = {'metabolites', ..};
%       me_options.id.reactions_kinetic_models   = {'actions', ..};
%     me_options.fba_constraints   (refer to the network model; are extended to augmented network below)
%
%Remarks:
% kinetic models in list 'kinetic models' are ordered by priority: highest priority first
% if the network is automatically extended, the FBA constraints have to be given for the 
% final extended network
%
% Outputs
%   network
%   network_CoHid
%   network_combined
%   res

me_options_default = me_default_options(length(kinetic_models));
me_options         = join_struct(me_options_default, me_options);
if isempty(me_options.id),
  me_options.id = me_options_default.id;
end

% -------------------------------------------------------------------------------
% Preparations 

% All compartments in kinetic models must also appear in network model
% and have the same sizes

display('  Checking the use of compartments');

if ~isfield(network,'compartments'),
  network.compartments      = {};
  network.compartment_sizes = [];
end

for it = 1:length(kinetic_models),
  if ~isfield(kinetic_models{it},'compartments'),
    kinetic_models{it}.compartments = {};
    kinetic_models{it}.compartment_sizes = [];
  end
  missing_compartments = setdiff(kinetic_models{it}.compartments,network.compartments);
  if length(missing_compartments),
    error(sprintf('Some compartments from model %d are missing in network model',it));
  end
  ll = label_names(kinetic_models{it}.compartments,network.compartments);
  if sum(kinetic_models{it}.compartment_sizes(find(ll)) ~= network.compartment_sizes(ll(find(ll)))),
    kinetic_models{it}.compartment_sizes(find(ll))
    network.compartment_sizes(ll(find(ll)))
    error('Compartment sizes do not match between models');
  end
end


% -------------------------------------------------------------------------------
% Phase 1: Map the kinetic models onto the network model

% Match elements and, if necessary, 
% add metabolites and reactions to the network model

[mapping_metabolites, mapping_reactions, covered_metabolites, covered_reactions, shared_metabolites, shared_reactions, network_aug, network_aug_CoHid,network_aug_CoSplit] = embedding_element_mapping(kinetic_models,network,me_options);

network_aug.external(label_names(me_options.set_external,network_aug.metabolites) ) = 1;
network_aug.external(label_names(me_options.set_internal,network_aug.metabolites) ) = 0;


% -------------------------------------------------------------------------------
% Phase 2: Determine consistent fluxes

% o determine concentrations and fluxes within the kinetic models
% o determine equilibrium concentrations within kinetic models
% o kinetic models are considered one after the other (corresponding to their
%   priority order: each model can predefine concentrations for the subsequent models!

[collect_v_all,collect_v_kinetic,c_stat,v_stat,collect_c,collect_v,collect_mu,mu] = model_embedding_consistent_fluxes(network_aug,kinetic_models,me_options,mapping_metabolites,mapping_reactions);

for it = 1:length(kinetic_models),
  kinetic_models{it}.s_init = c_stat{it};
end


% ---------------------------------------
% project fluxes on network model, where fluxes are constrained fluxes on the kinetic model parts

% scale network fluxes to make them comparable to kinetic fluxes
ind_v_kinetic = find(isfinite(collect_v_kinetic));
ratio         = median(collect_v_kinetic(ind_v_kinetic) ./ collect_v_all(ind_v_kinetic) );

v_mean = ratio * collect_v_all;
v_std  = guess_flux_std(v_mean);
v_sign = nan * v_mean;
v_fix  = nan * v_mean;

nm = length(network.actions);

if isstruct(me_options.fba_constraints),
  v_sign(1:nm) = me_options.fba_constraints.v_sign;
  v_fix(1:nm)  = me_options.fba_constraints.v_fix;
end

% impose fluxes from kinetic model(s) as fixed constraint
v_fix(isfinite(collect_v_kinetic)) = collect_v_kinetic(isfinite(collect_v_kinetic));

if ~isfield(network_aug,'reaction_names'),
  network_aug.reaction_names = network_aug.actions;
end

try 
  v_proj = project_fluxes(network_aug.N, find(network_aug.external), v_mean, v_std, v_sign, struct('method','euclidean'),v_fix);
catch
  warning('Steady state fluxes of kinetic models cannot be embedded into the larger network; I am now trying to find an approximate solution');
  v_mean(find(isfinite(v_fix))) = v_fix(find(isfinite(v_fix)));
  v_std  = guess_flux_std(v_mean);
  v_proj = project_fluxes(network_aug.N, find(network_aug.external), v_mean, v_std, v_sign, struct('method','euclidean'));
end

% ---------------------------------------
% Determine enyzme adjustments for the kinetic models
% for each kinetic model: v[stationary within network] = enzyme_adjustment * v[original]

for it = 1:length(kinetic_models);
  enzyme_adjustment{it} = v_proj(mapping_reactions{it}) ./ v_stat{it};
  if sum(isnan(enzyme_adjustment{it})),
    error('Vanishing fluxes in submodel');
  end
  if min(abs(enzyme_adjustment{it}))<10^-5,
    error('Difference in fluxes too extreme for enzyme adjustment');
  end  
end


% -----------------------------------------------------------------------
% collect concentrations and chemical potentials from kinetic models

c_fix  = nan * ones(size(network_aug.metabolites));

for it = length(kinetic_models):-1:1;
  c_fix(mapping_metabolites{it}) = c_stat{it};
end

% Determine predefined dmu within each kinetic model; 
% copy them into a vector dmu_fix for the entire network

dmu_fix = nan * ones(size(network_aug.actions));
for it = 1:length(kinetic_models);
  mmm = mu{it};
  mmm(isnan(mmm)) = 0;
  my_dmu = kinetic_models{it}.N'* mmm;
  my_dmu(find(abs(kinetic_models{it}.N') * isnan(mu{it}))) = nan;
  dmu_fix(mapping_reactions{it}) = my_dmu;
end

% Check if the flux distribution is thermo-feasible, and determine a feasible dmu vector

[eba_feasible, dmu] = eba_feasible_lp(v_proj, network_aug.N, dmu_fix,[],[],8); 

if ~eba_feasible, 
  warning('Flux distribution is infeasible for given chemical potentials.');
  display('Workaround: Inventing chemical potential differences');
  %% ALTERNATIVE: INVENT dmu .. THIS HAS TO BE FIXED;
  %% PROBLEMS: original models may be infeasible; 
  %% numeric calculation of stationary and equilibrium concentrations 
  %% is not accurate enough to obtain dmu that agree with rate signs 
  [eba_feasible, dmu] = eba_feasible_lp(v_proj, network_aug.N,[],100); 
  if ~eba_feasible,
    warning('Infeasible flux distribution');
    display('Next workaround: Inventing chemical potential differences');
    dmu = -sign(v_proj);
  end
end


% -------------------------------------------------------------------------------
% Phase 3: Insert standard rate laws and build the combined model

% Problem to be solved: predefine the c and mu values 
% of communicating metabolites in elasticity sampling 

% run elasticity sampling

[es_options, es_constraints] = es_default_options(network_aug);

es_constraints.v_fix         = v_proj; 
es_constraints.log_c_fix     = log(c_fix); 
es_constraints.dmu_fix       = dmu;
es_constraints.dmu_limit     = 10;      
es_constraints.dmu_limit_min = 10^-5;

result = es_reference_state(network_aug, es_options, es_constraints);

if me_options.verbose,
  display('  Check: Predefined and resulting concentrations');
  [c_fix,  result.c]
  display('  Check: Predefined and resulting chem. pot differences');
  [dmu,   -result.A]
end


% -------------------------------------------------------------------------------
% Outputs

network_combined                              = network_aug;
network_combined.kinetics                     = struct;
network_combined.kinetics.type                = 'embedded_kinetic_models';
network_combined.kinetics.kinetic_models      = kinetic_models;
network_combined.kinetics.kinetics_network    = result.kinetics;
network_combined.kinetics.mapping_metabolites = mapping_metabolites;
network_combined.kinetics.mapping_reactions   = mapping_reactions;
network_combined.kinetics.enzyme_adjustment   = enzyme_adjustment;
network_combined.s_init                       = result.c;

res.v_mean              = v_mean;
res.c_fix               = c_fix;
res.dmu_fix             = dmu_fix;
res.mu                  = result.mu;
res.A                   = result.A;
res.c_stat              = c_stat;
res.v_stat              = v_stat;
res.v_proj              = v_proj;
res.c_combined          = network_combined.kinetics.kinetics_network.c;
res.v_combined          = network_velocities(res.c_combined,network_combined);
res.collect_c           = collect_c;
res.collect_v           = collect_v;
res.collect_mu          = collect_mu;
res.covered_reactions   = covered_reactions;
res.covered_metabolites = covered_metabolites;
res.shared_metabolites  = shared_metabolites;
res.shared_reactions    = shared_reactions;
res.mapping_metabolites = mapping_metabolites;
res.mapping_reactions   = mapping_reactions;
