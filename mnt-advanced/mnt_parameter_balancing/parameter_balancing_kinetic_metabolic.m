function [r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std] = parameter_balancing_kinetic_metabolic(network, kinetic_data, options, v);

%  [r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std] = parameter_balancing_kinetic(network, kinetic_data, file_kinetic_data, options, v);
%
% Determine consistent kinetic parameter and metabolic state by parameter balancing
%
% This function resembles parameter_balancing_kinetic. Please see the documentation there.

  epsilon = 0.00001;
  
% ------------------------------------------------------------------------
% Initialise
  
eval(default('kinetic_data', '[]', 'options', 'struct'));

options = join_struct(parameter_balancing_default_options, options);

[nm,nr] = size(network.N);


% ----------------------------------------------------------------
% Lists of quantities to be considered
% (consider using 'parameter_balancing_quantities' instead)

basic_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant','concentration','concentration of enzyme'}';

model_quantities  = {'standard chemical potential','catalytic rate constant geometric mean', 'Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant', 'Michaelis constant product', 'concentration','reaction affinity','concentration of enzyme'}';

data_quantities   = {'standard chemical potential', 'standard Gibbs energy of reaction', 'Michaelis constant','activation constant', 'inhibitory constant', 'equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant', 'concentration','reaction affinity','concentration of enzyme'}';


% ------------------------------------------------------------------------
% Load and modify parameter priors

parameter_prior = parameter_balancing_prior([],options.parameter_prior_file);
parameter_prior = pb_parameter_prior_adjust(parameter_prior, options); 


% ----------------------------------------------------------------
% Load and modify kinetic data

if isstr(kinetic_data),
  %% Load data from file
  kinetic_data = data_integration_load_kinetic_data(data_quantities, [], network, kinetic_data, 0, 1);
elseif isempty(kinetic_data),  
  %% If necessary, create empty kinetic_data structure
  kinetic_data = data_integration_load_kinetic_data(data_quantities, [], network, [], 0, 1, options.reaction_column_name, options.compound_column_name);
end  

kinetic_data_orig = kinetic_data;
kinetic_data      = pb_kinetic_data_adjust(kinetic_data, parameter_prior, network, options);

if options.enforce_flux_directions, 
  display('Enforcing predefined flux directions');
  kinetic_data.A.lower(v>0) = epsilon;
  kinetic_data.A.upper(v<0) = -epsilon;
end

% ----------------------------------------------------------------
% Run parameter balancing

network.kinetics            = set_kinetics(network, 'cs');
task                        = parameter_balancing_task(network, kinetic_data, parameter_prior, model_quantities, basic_quantities);
res                         = parameter_balancing_calculation(task, parameter_prior, options);
[r,r_mean,r_std,r_orig,r_samples]  = parameter_balancing_output(res, kinetic_data_orig, options);


if options.postprocessing_enforce_fluxes,
   display(sprintf('The balanced kcat values and enzyme levels have been adjusted to the predefined fluxes. Assuming rate law of type "%s"',r.type));
   display('  (Note that weighting with standard deviations has not yet been implemented)');
   display('  (Adjusted enzyme levels may be outside the allowed range)');
   v_pred = network_velocities(r.c,network,r);
  
   ind_wrong_direction = find(v .* v_pred < 0);
   if ind_wrong_direction,
     network.actions(ind_wrong_direction)
     error('Flux directions do not match predefined fluxes in the above reactions'); 
   end
   v_pred(v==0) = nan; 
   v(v==0)      = nan; 
   %% figure(100);
   %% plot(abs(v),abs(v_pred),'.','MarkerSize',20); set(gca,'XScale','log','YScale','log');
   %% xlabel('Fluxes (predefined)');    xlabel('Fluxes (predicted from balanced steady state)');
   scaling_factors = sqrt(abs(v) ./ abs(v_pred))
   ind = isfinite(scaling_factors);
   r.KV(ind)    = r.KV(ind)    .* scaling_factors(ind);
   r.Kcatf(ind) = r.Kcatf(ind) .* scaling_factors(ind);
   r.Kcatr(ind) = r.Kcatr(ind) .* scaling_factors(ind);
   r.u(ind)     = r.u(ind)     .* scaling_factors(ind);
   r_std.KV(ind)    = r_std.KV(ind)    .* scaling_factors(ind);
   r_std.Kcatf(ind) = r_std.Kcatf(ind) .* scaling_factors(ind);
   r_std.Kcatr(ind) = r_std.Kcatr(ind) .* scaling_factors(ind);
   r_std.u(ind)     = r_std.u(ind)     .* scaling_factors(ind);
end
