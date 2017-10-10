function [r_median, r_mean, r_std, r_orig, r_samples]  = parameter_balancing_output(res,kinetic_data_orig,options)
  
% [r_median, r_mean, r_std,r_orig,r_samples]  = parameter_balancing_output(res,kinetic_data_orig,options)
%
% Auxiliary function

eval(default('kinetic_data_orig','[]', 'options','struct'));

% ------------------------------------------------
% r: "kinetics" data structure, balanced values

r_median.type = options.kinetics;
r_median.h    = ones(size(res.kinetics_posterior_median.Kcatf));

if isfield(res.kinetics_posterior_median,'mu0'), r_median.mu0   = res.kinetics_posterior_median.mu0  ; end
if isfield(res.kinetics_posterior_median,'KV'),  r_median.KV    = res.kinetics_posterior_median.KV   ; end
if isfield(res.kinetics_posterior_median,'Keq'), r_median.Keq   = res.kinetics_posterior_median.Keq  ; end
r_median.KM    =  res.kinetics_posterior_median.KM   ;
r_median.KA    =  res.kinetics_posterior_median.KA   ;
r_median.KI    =  res.kinetics_posterior_median.KI   ;
r_median.Kcatf =  res.kinetics_posterior_median.Kcatf;
r_median.Kcatr =  res.kinetics_posterior_median.Kcatr;
if isfield(res.kinetics_posterior_median,'c'),  r_median.c   = res.kinetics_posterior_median.c;  end
if isfield(res.kinetics_posterior_median,'u'),  r_median.u   = res.kinetics_posterior_median.u;  end

% ------------------------------------------------
% r_mean

if isfield(res.kinetics_posterior_mean,'mu0'), r_mean.mu0   = res.kinetics_posterior_mean.mu0  ; end
if isfield(res.kinetics_posterior_mean,'KV'),  r_mean.KV    = res.kinetics_posterior_mean.KV   ; end
if isfield(res.kinetics_posterior_mean,'Keq'), r_mean.Keq   = res.kinetics_posterior_mean.Keq  ; end
r_mean.KM    =  res.kinetics_posterior_mean.KM   ;
r_mean.KA    =  res.kinetics_posterior_mean.KA   ;
r_mean.KI    =  res.kinetics_posterior_mean.KI   ;
r_mean.Kcatf =  res.kinetics_posterior_mean.Kcatf;
r_mean.Kcatr =  res.kinetics_posterior_mean.Kcatr;
if isfield(res.kinetics_posterior_mean,'c'),  r_mean.c   = res.kinetics_posterior_mean.c;  end
if isfield(res.kinetics_posterior_mean,'u'),  r_mean.u   = res.kinetics_posterior_mean.u;  end

% ------------------------------------------------
% r_std: standard deviations

if isfield(res.kinetics_posterior_std,'mu0'),  r_std.mu0   = full(res.kinetics_posterior_std.mu0)  ; end
if isfield(res.kinetics_posterior_std,'KV'),   r_std.KV    = full(res.kinetics_posterior_std.KV)   ;  end
r_std.KM = res.kinetics_posterior_std.KM   ;
r_std.KA = res.kinetics_posterior_std.KA   ;
r_std.KI = res.kinetics_posterior_std.KI   ;
if isfield(res.kinetics_posterior_std,'Keq'),  r_std.Keq   = full(res.kinetics_posterior_std.Keq  ); end
r_std.Kcatf = full(res.kinetics_posterior_std.Kcatf);
r_std.Kcatr = full(res.kinetics_posterior_std.Kcatr);
if isfield(res.kinetics_posterior_std,'c'),   r_std.c     = res.kinetics_posterior_std.c; end
if isfield(res.kinetics_posterior_std,'u'),   r_std.u     = res.kinetics_posterior_std.u; end

% ------------------------------------------------
% r_orig: original values

if length(kinetic_data_orig),
  
  r_orig.mu0   = kinetic_data_orig.mu0.median  ;
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
