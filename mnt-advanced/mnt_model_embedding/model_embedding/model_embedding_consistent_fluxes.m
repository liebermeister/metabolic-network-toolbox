function [collect_v_all,collect_v_kinetic,c_stat,v_stat,collect_c,collect_v,collect_mu,mu] = model_embedding_consistent_fluxes(network,kinetic_models,me_options,mapping_metabolites,mapping_reactions)

% o determine concentrations and fluxes within the kinetic models
% o determine equilibrium concentrations within kinetic models
% o kinetic models are considered one after the other (corresponding to their
%   priority order: each model can predefine concentrations for the subsequent models!

% values resulting from each model's own concentrations
collect_c_own= nan * ones(length(network.metabolites),length(kinetic_models));
collect_v_own= nan * ones(length(network.actions),length(kinetic_models));

% values possibly changed by imposed concentrations
collect_c    = nan * ones(length(network.metabolites),length(kinetic_models));
collect_v    = nan * ones(length(network.actions),length(kinetic_models));
collect_c_eq = nan * ones(length(network.metabolites),length(kinetic_models));
collect_mu   = nan * ones(length(network.metabolites),length(kinetic_models));

% collects all values fixed so far
collect_c_all    = nan * ones(length(network.metabolites),1);
collect_c_eq_all = nan * ones(length(network.metabolites),1);
collect_v_all    = nan * ones(length(network.actions),1);
collect_v_kinetic= nan * ones(length(network.actions),1);

% insert fluxes for network (if given)
if isfield(me_options,'v_network'),
  collect_v_all(1:length(me_options.v_network)) = me_options.v_network;
end

for it = 1:length(kinetic_models),

  %% concentrations and rates as given by kinetic models individually
  nn              = kinetic_models{it};
  c_stat_own{it}  = nn.s_init;
  v_stat_own{it}  = network_velocities(c_stat_own{it},nn);
  collect_c_own(mapping_metabolites{it},it) = c_stat_own{it};
  collect_v_own(mapping_reactions{it},it)     = v_stat_own{it};

  %% ... but some concentrations get fixed by preceding models!
  ind_predefined_c = find(isfinite(collect_c_all(mapping_metabolites{it})));
  predefined_c     = collect_c_all(mapping_metabolites{it}(ind_predefined_c));
  predefined_c_eq  = collect_c_eq_all(mapping_metabolites{it}(ind_predefined_c));
  nn.external(ind_predefined_c) = 1;
  c_stat{it} = c_stat_own{it};
  c_stat{it}(ind_predefined_c) = predefined_c;
  v_stat{it} = network_velocities(c_stat{it}, kinetic_models{it});

  if me_options.make_kinetic_models_stationary,
    c_stat{it}  = network_steady_state(nn,c_stat{it},10000);
    v_stat{it}  = network_velocities(c_stat{it}, kinetic_models{it});
    v_stat{it}  = es_make_fluxes_stationary(kinetic_models{it},  v_stat{it});
    v_stat{it}  = es_make_fluxes_stationary(kinetic_models{it},  v_stat{it});
  end

  %% compute equilibrium concentrations (only possible for reversible rate laws!!!
  nnn          = nn;
  nnn.external = 0 * nnn.external;
  c_eq{it}     = c_stat{it};
  c_eq{it}(ind_predefined_c) = predefined_c_eq;
  c_eq{it}    = network_steady_state(nnn,  c_eq{it},10);%000);

  %% for safety reasons (in case of irreversible rate laws)
  c_eq{it} = real(c_eq{it});
  c_eq{it}(c_eq{it}<10^-15) = 10^-15;
  
  %% compare rate signs and driving forces
  %[ sign(v_stat{it}) .* sign(-kinetic_models{it}.N'*log(c_stat{it}./c_eq{it}))]
  mu{it}   = RT * log(c_stat{it}./c_eq{it});

  collect_c(mapping_metabolites{it},it)    = c_stat{it};
  collect_v(mapping_reactions{it},it)      = v_stat{it};
  collect_c_eq(mapping_metabolites{it},it) = c_eq{it};
  collect_mu(mapping_metabolites{it},it)   = mu{it};
  
  ind_undet = find(isnan(collect_c_all));
  collect_c_all(ind_undet) = collect_c(ind_undet,it);
  collect_c_eq_all(ind_undet) = collect_c_eq(ind_undet,it);
  ind_undet = find(isnan(collect_v_all));
  collect_v_all(ind_undet)     = collect_v(ind_undet,it);
  
  ind_det = find(isfinite(collect_v(:,it)));
  collect_v_kinetic(ind_det) = collect_v(ind_det,it);

end

% ---------------------------------------
% if there are overlaps between kinetic models and concentrations or equilibrium concentrations differ -> error

for it = 1:length(kinetic_models);
  if find([c_eq{it}==0]+[isfinite(c_eq{it})==0]), warning('Irreversible reaction encountered'); end
end

for it = 1:length(network.metabolites),
  dum = isfinite(collect_c_eq(it,:));
  if length(unique(collect_c_eq(it,dum))) > 1, 
    warning(sprintf('Inconsistent equilibrium concentration of metabolite %s',network.metabolites{it})); 
    % collect_c_eq
  end
end
