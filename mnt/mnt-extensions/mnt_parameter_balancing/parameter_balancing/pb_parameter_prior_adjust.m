function parameter_prior = pb_parameter_prior_adjust(parameter_prior, options)

% parameter_prior = pb_parameter_prior_adjust(parameter_prior, options)
%
% Replace same values in prior table according to user options:
% 
% options.kcat_prior_median
% options.kcat_prior_log10_std
% options.kcat_lower
% options.kcatr_lower
% options.KM_lower
% options.keq_upper
% options.kcat_upper
% options.conc_min
% options.conc_max

if length(options.kcat_prior_median),
  display(sprintf('  Adjusting the prior: Setting kcat prior median value to %f', options.kcat_prior_median));
  ll = label_names({'catalytic rate constant geometric mean', 'substrate catalytic rate constant','product catalytic rate constant'},  parameter_prior.QuantityType);
  parameter_prior.PriorMedian(ll) = repmat({num2str(options.kcat_prior_median)},3,1);
end

if length(options.kcat_prior_log10_std),
  display(sprintf('  Adjusting the prior: Setting kcat prior log10 std dev value to %f', options.kcat_prior_log10_std));
  ll = label_names({'catalytic rate constant geometric mean', 'substrate catalytic rate constant','product catalytic rate constant'},  parameter_prior.QuantityType);
  parameter_prior.PriorStd(ll) = repmat({num2str(options.kcat_prior_log10_std)},3,1);
end

if length(options.kcat_lower),
  display(sprintf('  Adjusting the prior: Setting kcat lower bound to %f', options.kcat_lower));
  ll = label_names({'catalytic rate constant geometric mean', 'substrate catalytic rate constant','product catalytic rate constant'},  parameter_prior.QuantityType);
  parameter_prior.LowerBound(ll) = repmat({num2str(options.kcat_lower)},3,1);
end

if length(options.kcatr_lower),
  display(sprintf('  Adjusting the prior: Setting kcat reverse lower bound to %f', options.kcatr_lower));
  ll = label_names({'product catalytic rate constant'},  parameter_prior.QuantityType);
  parameter_prior.LowerBound(ll) = repmat({num2str(options.kcatr_lower)},1,1);
end

if length(options.KM_lower),
  display(sprintf('  Adjusting the prior: Setting KM lower bound to %2.3f', options.KM_lower));
  ll = label_names({'Michaelis constant'},  parameter_prior.QuantityType);
  parameter_prior.LowerBound(ll) = {num2str(options.KM_lower)};
end

if length(options.KM_lower),
  display(sprintf('  Adjusting the prior: Setting KM upper bound to %2.3f', options.KM_upper));
  ll = label_names({'Michaelis constant'},  parameter_prior.QuantityType);
  parameter_prior.UpperBound(ll) = {num2str(options.KM_upper)};
end

if length(options.Keq_upper),
  display(sprintf('  Adjusting the prior: Setting Keq upper bound to %2.3f', options.Keq_upper));
  ll = label_names({'equilibrium constant'},  parameter_prior.QuantityType);
  parameter_prior.UpperBound(ll) = {num2str(options.Keq_upper)};
end

if length(options.kcat_upper),
  display(sprintf('  Adjusting the prior: Setting kcat upper bound to %2.3f', options.kcat_upper));
  ll = label_names({'catalytic rate constant geometric mean', 'substrate catalytic rate constant','product catalytic rate constant'},  parameter_prior.QuantityType);
  parameter_prior.UpperBound(ll) = repmat({num2str(options.kcat_upper)},3,1);
end

ll = label_names({'concentration'},  parameter_prior.QuantityType);

if options.conc_min ~= eval(parameter_prior.LowerBound{ll}),
  display(sprintf('Setting default minimum concentration to %f',options.conc_min)); 
  parameter_prior.LowerBound{ll} = num2str(options.conc_min); 
end

if options.conc_max ~= eval(parameter_prior.UpperBound{ll}), 
  display(sprintf('Setting default maximum concentration to %f',options.conc_max)); 
  parameter_prior.UpperBound{ll} = num2str(options.conc_max); 
end
