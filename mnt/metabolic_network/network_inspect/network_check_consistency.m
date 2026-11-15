function network_check_consistency(network)

% network_check_consistency(network)
%
% Checks whether the network admits feasible fluxes. Biomass reactions are ignored: the reactions in
% network.biomass_reactions if this field exists, otherwise a reaction called 'Biomass production'.

[options,constraints]   = flux_constraints_default_options(network);
ind_biomass_production  = find(ismember(network.actions, network_biomass_reactions(network, {'Biomass production'})));
constraints.ind_ignore  = ind_biomass_production;
J                       = sample_feasible_v(network.N,find(network.external),constraints,options);

if length(J)
  display('There are feasible fluxes');
else
  display('There are no feasible fluxes');
end


% check existence of equilibrium reactions
ind_eq_reactions = network_equilibrium_reactions(network,1);

if ~length(ind_eq_reactions)
  display('There are no structurally determined equilibrium reactions');
end
