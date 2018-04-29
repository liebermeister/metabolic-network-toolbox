function [r_mode, r_mean, r_std, r_geom_mean, r_geom_std, r_orig, r_samples]  = parameter_balancing_output(res, kinetic_data_orig, options)
  
% [r_mode, r_mean, r_std, r_geom_mean, r_geom_std, r_orig,r_samples]  = parameter_balancing_output(res,kinetic_data_orig,options)
%
% Auxiliary function

eval(default('kinetic_data_orig','[]', 'options','struct'));

% ------------------------------------------------
% r: "kinetics" data structure, balanced values

r_mode.type  = options.kinetics;
r_mode.h     = ones(size(res.kinetics.posterior_mode.Kcatf));
r_mode.KM    =  res.kinetics.posterior_mode.KM   ;
r_mode.KA    =  res.kinetics.posterior_mode.KA   ;
r_mode.KI    =  res.kinetics.posterior_mode.KI   ;
r_mode.Kcatf =  res.kinetics.posterior_mode.Kcatf;
r_mode.Kcatr =  res.kinetics.posterior_mode.Kcatr;

if isfield(res.kinetics.posterior_mode,'mu0'), r_mode.mu0   = res.kinetics.posterior_mode.mu0  ; end
if isfield(res.kinetics.posterior_mode,'mu'),  r_mode.mu    = res.kinetics.posterior_mode.mu   ; end
if isfield(res.kinetics.posterior_mode,'KV'),  r_mode.KV    = res.kinetics.posterior_mode.KV   ; end
if isfield(res.kinetics.posterior_mode,'Keq'), r_mode.Keq   = res.kinetics.posterior_mode.Keq  ; end
if isfield(res.kinetics.posterior_mode,'c'),   r_mode.c     = res.kinetics.posterior_mode.c;  end
if isfield(res.kinetics.posterior_mode,'u'),   r_mode.u     = res.kinetics.posterior_mode.u;  end
if isfield(res.kinetics.posterior_mode,'A'),   r_mode.A     = res.kinetics.posterior_mode.A;  end

% ------------------------------------------------
% r_mean

r_mean.KM    =  res.kinetics.unconstrained_posterior_mean.KM   ;
r_mean.KA    =  res.kinetics.unconstrained_posterior_mean.KA   ;
r_mean.KI    =  res.kinetics.unconstrained_posterior_mean.KI   ;
r_mean.Kcatf =  res.kinetics.unconstrained_posterior_mean.Kcatf;
r_mean.Kcatr =  res.kinetics.unconstrained_posterior_mean.Kcatr;
if isfield(res.kinetics.unconstrained_posterior_mean,'mu0'), r_mean.mu0   = res.kinetics.unconstrained_posterior_mean.mu0  ; end
if isfield(res.kinetics.unconstrained_posterior_mean,'mu'),  r_mean.mu    = res.kinetics.unconstrained_posterior_mean.mu  ; end
if isfield(res.kinetics.unconstrained_posterior_mean,'KV'),  r_mean.KV    = res.kinetics.unconstrained_posterior_mean.KV   ; end
if isfield(res.kinetics.unconstrained_posterior_mean,'Keq'), r_mean.Keq   = res.kinetics.unconstrained_posterior_mean.Keq  ; end
if isfield(res.kinetics.unconstrained_posterior_mean,'c'),   r_mean.c     = res.kinetics.unconstrained_posterior_mean.c;  end
if isfield(res.kinetics.unconstrained_posterior_mean,'u'),   r_mean.u     = res.kinetics.unconstrained_posterior_mean.u;  end
if isfield(res.kinetics.unconstrained_posterior_mean,'A'),   r_mean.A     = res.kinetics.unconstrained_posterior_mean.A;  end

% ------------------------------------------------
% r_std: standard deviations

r_std.KM    = res.kinetics.unconstrained_posterior_std.KM   ;
r_std.KA    = res.kinetics.unconstrained_posterior_std.KA   ;
r_std.KI    = res.kinetics.unconstrained_posterior_std.KI   ;
r_std.Kcatf = full(res.kinetics.unconstrained_posterior_std.Kcatf);
r_std.Kcatr = full(res.kinetics.unconstrained_posterior_std.Kcatr);
if isfield(res.kinetics.unconstrained_posterior_std,'mu0'),  r_std.mu0   = full(res.kinetics.unconstrained_posterior_std.mu0)  ; end
if isfield(res.kinetics.unconstrained_posterior_std,'mu'),   r_std.mu    = full(res.kinetics.unconstrained_posterior_std.mu)   ; end
if isfield(res.kinetics.unconstrained_posterior_std,'KV'),   r_std.KV    = full(res.kinetics.unconstrained_posterior_std.KV)   ;  end
if isfield(res.kinetics.unconstrained_posterior_std,'Keq'),  r_std.Keq   = full(res.kinetics.unconstrained_posterior_std.Keq  ); end
if isfield(res.kinetics.unconstrained_posterior_std,'c'),   r_std.c     = res.kinetics.unconstrained_posterior_std.c; end
if isfield(res.kinetics.unconstrained_posterior_std,'u'),   r_std.u     = res.kinetics.unconstrained_posterior_std.u; end
if isfield(res.kinetics.unconstrained_posterior_std,'A'),   r_std.A     = res.kinetics.unconstrained_posterior_std.A; end

% ------------------------------------------------
% r_geom_mean: geometric mean

r_geom_mean.KM    = res.kinetics.unconstrained_posterior_median.KM   ;
r_geom_mean.KA    = res.kinetics.unconstrained_posterior_median.KA   ;
r_geom_mean.KI    = res.kinetics.unconstrained_posterior_median.KI   ;
r_geom_mean.Kcatf = full(res.kinetics.unconstrained_posterior_median.Kcatf);
r_geom_mean.Kcatr = full(res.kinetics.unconstrained_posterior_median.Kcatr);
if isfield(res.kinetics.unconstrained_posterior_median,'mu0'),  r_geom_mean.mu0   = full(res.kinetics.unconstrained_posterior_median.mu0)  ; end
if isfield(res.kinetics.unconstrained_posterior_median,'mu'),   r_geom_mean.mu    = full(res.kinetics.unconstrained_posterior_median.mu)   ; end
if isfield(res.kinetics.unconstrained_posterior_median,'KV'),   r_geom_mean.KV    = full(res.kinetics.unconstrained_posterior_median.KV)   ;  end
if isfield(res.kinetics.unconstrained_posterior_median,'Keq'),  r_geom_mean.Keq   = full(res.kinetics.unconstrained_posterior_median.Keq  ); end
if isfield(res.kinetics.unconstrained_posterior_median,'c'),    r_geom_mean.c     = res.kinetics.unconstrained_posterior_median.c; end
if isfield(res.kinetics.unconstrained_posterior_median,'u'),    r_geom_mean.u     = res.kinetics.unconstrained_posterior_median.u; end
if isfield(res.kinetics.unconstrained_posterior_median,'A'),    r_geom_mean.A     = res.kinetics.unconstrained_posterior_median.A; end

% ------------------------------------------------
% r_std: geometric standard deviations

r_geom_std.KM    = res.kinetics.unconstrained_posterior_geom_std.KM   ;
r_geom_std.KA    = res.kinetics.unconstrained_posterior_geom_std.KA   ;
r_geom_std.KI    = res.kinetics.unconstrained_posterior_geom_std.KI   ;
r_geom_std.Kcatf = full(res.kinetics.unconstrained_posterior_geom_std.Kcatf);
r_geom_std.Kcatr = full(res.kinetics.unconstrained_posterior_geom_std.Kcatr);
if isfield(res.kinetics.unconstrained_posterior_geom_std,'mu0'),  r_geom_std.mu0   = full(res.kinetics.unconstrained_posterior_geom_std.mu0)  ; end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'mu'),   r_geom_std.mu    = full(res.kinetics.unconstrained_posterior_geom_std.mu)   ; end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'KV'),   r_geom_std.KV    = full(res.kinetics.unconstrained_posterior_geom_std.KV)   ;  end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'Keq'),  r_geom_std.Keq   = full(res.kinetics.unconstrained_posterior_geom_std.Keq  ); end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'c'),   r_geom_std.c     = res.kinetics.unconstrained_posterior_geom_std.c; end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'u'),   r_geom_std.u     = res.kinetics.unconstrained_posterior_geom_std.u; end
if isfield(res.kinetics.unconstrained_posterior_geom_std,'A'),   r_geom_std.A     = res.kinetics.unconstrained_posterior_geom_std.A; end


% ------------------------------------------------
% r_orig: original values

if length(kinetic_data_orig),
  
  r_orig.mu0   = kinetic_data_orig.mu0.median  ;
  if isfield(kinetic_data_orig,'mu'),
    r_orig.mu    = kinetic_data_orig.mu.median  ;
  end
  r_orig.KM    = kinetic_data_orig.KM.median   ;
  r_orig.KA    = kinetic_data_orig.KA.median   ;
  r_orig.KI    = kinetic_data_orig.KI.median   ;
  r_orig.Keq   = kinetic_data_orig.Keq.median  ;
  if isfield(kinetic_data_orig,'c'),
    r_orig.c     = kinetic_data_orig.c.median   ;
  end
  if isfield(kinetic_data_orig,'u'),
    r_orig.u     = kinetic_data_orig.u.median  ;
  end
  
  
  if isfield(options,'Keq_given'),
    if length(options.Keq_given),
      ind_finite = find(isfinite(options.Keq_given));
      r_orig.Keq(ind_finite) = options.Keq_given(ind_finite);
    end
  end
  
  if isfield(r_orig,'dmu0'), 
    r_orig.Keq(isnan(r_orig.Keq)) = exp(-1/RT * r_orig.dmu0(isnan(r_orig.Keq)));
  end


  r_orig.Kcatf = kinetic_data_orig.Kcatf.median;
  r_orig.Kcatr = kinetic_data_orig.Kcatr.median;
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
