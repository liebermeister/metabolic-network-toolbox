function [r_mode, r_mean, r_std, r_geom_mean, r_geom_std, r_orig, r_samples, r_flux_adjusted, r_flux_adjusted_std]  = parameter_balancing_output(res, kinetic_data_orig, options)
  
% [r_mode, r_mean, r_std, r_geom_mean, r_geom_std, r_orig, r_samples, r_flux_adjusted, r_flux_adjusted_std]  = parameter_balancing_output(res,kinetic_data_orig,options)
%
% Convert parameter sets, previously obtained by parameter balancing, into the format 
% used to describe kinetic parameters in the Metabolic Network Toolbox.
%
% If option "adjust_to_fluxes" is set, the resulting parameter set (constrained posterior mode)
% is additionally adjusted to yield the predefined flux distribution (vector given in options.v)
  
eval(default('kinetic_data_orig','[]', 'options','struct'));

% ------------------------------------------------
% r: "kinetics" data structure, balanced values

r_mode.type  = options.kinetics;
r_mode.h     = myround(ones(size(res.kinetics.posterior_mode.Kcatf)));
r_mode.KM    = myround(res.kinetics.posterior_mode.KM   );
r_mode.KA    = myround(res.kinetics.posterior_mode.KA   );
r_mode.KI    = myround(res.kinetics.posterior_mode.KI   );
r_mode.Kcatf = myround(res.kinetics.posterior_mode.Kcatf);
r_mode.Kcatr = myround(res.kinetics.posterior_mode.Kcatr);

if isfield(res.kinetics.posterior_mode,'mu0'), r_mode.mu0   = myround(res.kinetics.posterior_mode.mu0); end
if isfield(res.kinetics.posterior_mode,'mu'),  r_mode.mu    = myround(res.kinetics.posterior_mode.mu );  end
if isfield(res.kinetics.posterior_mode,'KV'),  r_mode.KV    = myround(res.kinetics.posterior_mode.KV );  end
if isfield(res.kinetics.posterior_mode,'Keq'), r_mode.Keq   = myround(res.kinetics.posterior_mode.Keq); end
if isfield(res.kinetics.posterior_mode,'c'),   r_mode.c     = myround(res.kinetics.posterior_mode.c  );   end
if isfield(res.kinetics.posterior_mode,'u'),   r_mode.u     = myround(res.kinetics.posterior_mode.u  );   end
if isfield(res.kinetics.posterior_mode,'A'),   r_mode.A     = myround(res.kinetics.posterior_mode.A  );   end

% ------------------------------------------------
% r_mean

r_mean.KM    =  myround(res.kinetics.unconstrained_posterior_mean.KM   );
r_mean.KA    =  myround(res.kinetics.unconstrained_posterior_mean.KA   );
r_mean.KI    =  myround(res.kinetics.unconstrained_posterior_mean.KI   );
r_mean.Kcatf =  myround(res.kinetics.unconstrained_posterior_mean.Kcatf);
r_mean.Kcatr =  myround(res.kinetics.unconstrained_posterior_mean.Kcatr);

if isfield(res.kinetics.unconstrained_posterior_mean,'mu0'), r_mean.mu0   = myround(res.kinetics.unconstrained_posterior_mean.mu0)      ;       end
if isfield(res.kinetics.unconstrained_posterior_mean,'mu'),  r_mean.mu    = myround(res.kinetics.unconstrained_posterior_mean.mu    );  end
if isfield(res.kinetics.unconstrained_posterior_mean,'KV'),  r_mean.KV    = myround(res.kinetics.unconstrained_posterior_mean.KV    );  end
if isfield(res.kinetics.unconstrained_posterior_mean,'Keq'), r_mean.Keq   = myround(res.kinetics.unconstrained_posterior_mean.Keq   ); end
if isfield(res.kinetics.unconstrained_posterior_mean,'c'),   r_mean.c     = myround(res.kinetics.unconstrained_posterior_mean.c     );   end
if isfield(res.kinetics.unconstrained_posterior_mean,'u'),   r_mean.u     = myround(res.kinetics.unconstrained_posterior_mean.u     );   end
if isfield(res.kinetics.unconstrained_posterior_mean,'A'),   r_mean.A     = myround(res.kinetics.unconstrained_posterior_mean.A     );   end

% ------------------------------------------------
% r_std: standard deviations

r_std.KM    = myround(res.kinetics.unconstrained_posterior_std.KM         );
r_std.KA    = myround(res.kinetics.unconstrained_posterior_std.KA         );
r_std.KI    = myround(res.kinetics.unconstrained_posterior_std.KI         );
r_std.Kcatf = myround(full(res.kinetics.unconstrained_posterior_std.Kcatf));
r_std.Kcatr = myround(full(res.kinetics.unconstrained_posterior_std.Kcatr));

if isfield(res.kinetics.unconstrained_posterior_std,'mu0'),  r_std.mu0   = myround(full(res.kinetics.unconstrained_posterior_std.mu0) ) ; end
if isfield(res.kinetics.unconstrained_posterior_std,'mu'),   r_std.mu    = myround(full(res.kinetics.unconstrained_posterior_std.mu)  ) ; end
if isfield(res.kinetics.unconstrained_posterior_std,'KV'),   r_std.KV    = myround(full(res.kinetics.unconstrained_posterior_std.KV)  ) ;  end
if isfield(res.kinetics.unconstrained_posterior_std,'Keq'),  r_std.Keq   = myround(full(res.kinetics.unconstrained_posterior_std.Keq) ); end
if isfield(res.kinetics.unconstrained_posterior_std,'c'),   r_std.c      = myround(res.kinetics.unconstrained_posterior_std.c         ); end
if isfield(res.kinetics.unconstrained_posterior_std,'u'),   r_std.u      = myround(res.kinetics.unconstrained_posterior_std.u         ); end
if isfield(res.kinetics.unconstrained_posterior_std,'A'),   r_std.A      = myround(res.kinetics.unconstrained_posterior_std.A         ); end

% ------------------------------------------------
% r_geom_mean: geometric mean

r_geom_mean.KM    = myround(res.kinetics.unconstrained_posterior_median.KM         )   ;
r_geom_mean.KA    = myround(res.kinetics.unconstrained_posterior_median.KA         )   ;
r_geom_mean.KI    = myround(res.kinetics.unconstrained_posterior_median.KI         )   ;
r_geom_mean.Kcatf = myround(full(res.kinetics.unconstrained_posterior_median.Kcatf));
r_geom_mean.Kcatr = myround(full(res.kinetics.unconstrained_posterior_median.Kcatr));

if isfield(res.kinetics.unconstrained_posterior_median,'mu0'),  r_geom_mean.mu0   = myround(full(res.kinetics.unconstrained_posterior_median.mu0))  ; end
if isfield(res.kinetics.unconstrained_posterior_median,'mu'),   r_geom_mean.mu    = myround(full(res.kinetics.unconstrained_posterior_median.mu) )  ; end
if isfield(res.kinetics.unconstrained_posterior_median,'KV'),   r_geom_mean.KV    = myround(full(res.kinetics.unconstrained_posterior_median.KV) )  ;  end
if isfield(res.kinetics.unconstrained_posterior_median,'Keq'),  r_geom_mean.Keq   = myround(full(res.kinetics.unconstrained_posterior_median.Keq)); end
if isfield(res.kinetics.unconstrained_posterior_median,'c'),    r_geom_mean.c     = myround(res.kinetics.unconstrained_posterior_median.c        ); end
if isfield(res.kinetics.unconstrained_posterior_median,'u'),    r_geom_mean.u     = myround(res.kinetics.unconstrained_posterior_median.u        ); end
if isfield(res.kinetics.unconstrained_posterior_median,'A'),    r_geom_mean.A     = myround(res.kinetics.unconstrained_posterior_median.A        ); end

% ------------------------------------------------
% r_std: geometric standard deviations

r_geom_std.KM    = myround(res.kinetics.unconstrained_posterior_geom_std.KM         ) ;
r_geom_std.KA    = myround(res.kinetics.unconstrained_posterior_geom_std.KA         ) ;
r_geom_std.KI    = myround(res.kinetics.unconstrained_posterior_geom_std.KI         ) ;
r_geom_std.Kcatf = myround(full(res.kinetics.unconstrained_posterior_geom_std.Kcatf));
r_geom_std.Kcatr = myround(full(res.kinetics.unconstrained_posterior_geom_std.Kcatr));

if isfield(res.kinetics.unconstrained_posterior_geom_std,'mu0'),  r_geom_std.mu0   = myround(full(res.kinetics.unconstrained_posterior_geom_std.mu0))  ; end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'mu'),   r_geom_std.mu    = myround(full(res.kinetics.unconstrained_posterior_geom_std.mu) )  ; end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'KV'),   r_geom_std.KV    = myround(full(res.kinetics.unconstrained_posterior_geom_std.KV) )  ;  end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'Keq'),  r_geom_std.Keq   = myround(full(res.kinetics.unconstrained_posterior_geom_std.Keq)); end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'c'),    r_geom_std.c     = myround(res.kinetics.unconstrained_posterior_geom_std.c        ); end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'u'),    r_geom_std.u     = myround(res.kinetics.unconstrained_posterior_geom_std.u        ); end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'A'),    r_geom_std.A     = myround(res.kinetics.unconstrained_posterior_geom_std.A        ); end


% ------------------------------------------------
% r_orig: original values

if length(kinetic_data_orig),
  
  %% if isfield(kinetic_data_orig,'mu'),
  %%   r_orig.mu0   = kinetic_data_orig.mu0.median  ;
  %%   r_orig.mu    = kinetic_data_orig.mu.median  ;
  %% end
  r_orig.KM    = myround(kinetic_data_orig.KM.median   );
  r_orig.KA    = myround(kinetic_data_orig.KA.median   );
  r_orig.KI    = myround(kinetic_data_orig.KI.median   );
  r_orig.Keq   = myround(kinetic_data_orig.Keq.median  );
  
  if isfield(kinetic_data_orig,'c'),    r_orig.c     = myround(kinetic_data_orig.c.median  ) ;   end
  if isfield(kinetic_data_orig,'u'),    r_orig.u     = myround(kinetic_data_orig.u.median  );  end
  if isfield(kinetic_data_orig,'KV'),   r_orig.KV    = myround(kinetic_data_orig.KV.median )  ;   end
  
  if isfield(options,'Keq_given'),
    if length(options.Keq_given),
      ind_finite = find(isfinite(options.Keq_given));
      r_orig.Keq(ind_finite) = myround(my_options.Keq_given(ind_finite));
    end
  end
  
  if isfield(r_orig,'dmu0'), 
    r_orig.Keq(isnan(r_orig.Keq)) = myround(exp(-1/RT * r_orig.dmu0(isnan(r_orig.Keq))));
  end

  r_orig.Kcatf = myround(kinetic_data_orig.Kcatf.median);
  r_orig.Kcatr = myround(kinetic_data_orig.Kcatr.median);
end


% ------------------------------------------------
% r_samples: sampled values

if isfield(options,'n_samples'),
  if options.n_samples,
    r_samples = res.kinetics_posterior_samples;
  else
    r_samples = [];
  end
end

% ----------------------------------------------------------------
% Adjust solution (in struct 'r_mode') to predefined fluxes

if options.adjust_to_fluxes,
  
  v = options.v;
  display(sprintf('  (The balanced kcat values and enzyme levels have been adjusted to the predefined fluxes.)\n  (Assuming rate law of type "%s")',r.type));
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
   scaling_factors = sqrt(abs(v) ./ abs(v_pred));
   ind = isfinite(scaling_factors);
   r_flux_adjusted = r;
   r_flux_adjusted_std = r_std;
   r_flux_adjusted.KV(ind)    = r_flux_adjusted.KV(ind)    .* scaling_factors(ind);
   r_flux_adjusted.Kcatf(ind) = r_flux_adjusted.Kcatf(ind) .* scaling_factors(ind);
   r_flux_adjusted.Kcatr(ind) = r_flux_adjusted.Kcatr(ind) .* scaling_factors(ind);
   r_flux_adjusted.u(ind)     = r_flux_adjusted.u(ind)     .* scaling_factors(ind);
   r_flux_adjusted_std.KV(ind)    = r_flux_adjusted_std.KV(ind)    .* scaling_factors(ind);
   r_flux_adjusted_std.Kcatf(ind) = r_flux_adjusted_std.Kcatf(ind) .* scaling_factors(ind);
   r_flux_adjusted_std.Kcatr(ind) = r_flux_adjusted_std.Kcatr(ind) .* scaling_factors(ind);
   r_flux_adjusted_std.u(ind)     = r_flux_adjusted_std.u(ind)     .* scaling_factors(ind);
end

% -------------------------------------------------

function myround = myround(x)

myround = round(x,4,'significant');
