function kinetic_data = kinetic_data_complete(kinetic_data,parameter_prior,insert_pseudo_values,network,verbose)

% kinetic_data = kinetic_data_complete(kinetic_data, parameter_prior, insert_pseudo_values,network);
%
% Complete missing entries in 'kinetic_data' data structure:
%  - complete missing logarithmic data values based on known non-logarithmic ones (and vice versa)
%  - if desired, replace unknown data values by pseudo values
%
% Prior values of equilibrium constants, standard chemical potentials and standard reaction Gibbs free
% energies are converted from the standard concentration of the prior table to 1 mM (the standard
% concentration used in kinetic_data); if the prior table does not use 1 mM, 'network' is required

eval(default('insert_pseudo_values','0','verbose','1','network','[]'));

if insert_pseudo_values
  if verbose
    display('  Replacing missing data values by pseudo values');
  end
end

fn = fieldnames(kinetic_data);

for it = 1:length(fn)

  clear my_lower my_upper

  ind_lower  = find(~isfinite(kinetic_data.(fn{it}).lower));
  ind_upper  = find(~isfinite(kinetic_data.(fn{it}).upper));
  my_scaling =  parameter_prior.MathematicalType{parameter_prior.symbol_index.(fn{it})};

  % shift from the standard concentration of the prior table to 1 mM (zero for most quantities)
  my_shift   = prior_standard_concentration_shift(parameter_prior, fn{it}, network);
  if isscalar(my_shift)
    my_shift = my_shift * ones(size(kinetic_data.(fn{it}).lower));
  end

  switch my_scaling

    case 'Multiplicative'

      % complete logarithmic bounds from given non-logarithmic bounds
      kinetic_data.(fn{it}) = complete_log_bounds(kinetic_data.(fn{it}));

      my_mean_ln   = log(eval(parameter_prior.PriorMedian{parameter_prior.symbol_index.(fn{it})}));
      my_std_ln    = log(10) * eval(parameter_prior.PriorStd{parameter_prior.symbol_index.(fn{it})});
      [my_mean, my_std] = lognormal_log2normal(my_mean_ln,my_std_ln);

      my_lower     = eval(parameter_prior.LowerBound{parameter_prior.symbol_index.(fn{it})});
      my_upper     = eval(parameter_prior.UpperBound{parameter_prior.symbol_index.(fn{it})});
      my_lower_ln  = log(my_lower);
      my_upper_ln  = log(my_upper);

      kinetic_data.(fn{it}).lower(ind_lower)    = my_lower * exp(my_shift(ind_lower));
      kinetic_data.(fn{it}).upper(ind_upper)    = my_upper * exp(my_shift(ind_upper));
      kinetic_data.(fn{it}).lower_ln(ind_lower) = my_lower_ln + my_shift(ind_lower);
      kinetic_data.(fn{it}).upper_ln(ind_upper) = my_upper_ln + my_shift(ind_upper);

    case 'Additive'

      my_mean    = eval(parameter_prior.PriorMedian{parameter_prior.symbol_index.(fn{it})});
      my_std     = eval(parameter_prior.PriorStd{parameter_prior.symbol_index.(fn{it})});

      my_lower     = eval(parameter_prior.LowerBound{parameter_prior.symbol_index.(fn{it})});
      my_upper     = eval(parameter_prior.UpperBound{parameter_prior.symbol_index.(fn{it})});
      kinetic_data.(fn{it}).lower(ind_lower) = my_lower + my_shift(ind_lower);
      kinetic_data.(fn{it}).upper(ind_upper) = my_upper + my_shift(ind_upper);
  end

  if insert_pseudo_values

    if strcmp('Derived',parameter_prior.Dependence{parameter_prior.symbol_index.(fn{it})})
      ind_mean  = find(~isfinite(kinetic_data.(fn{it}).mean));
      kinetic_data.(fn{it}).std(ind_mean)  = my_std;

      switch my_scaling
        case 'Additive'
          kinetic_data.(fn{it}).mean(ind_mean) = my_mean + my_shift(ind_mean);
          kinetic_data.(fn{it}).median = kinetic_data.(fn{it}).mean;
        case 'Multiplicative'
          % (mean and std of a log-normal distribution scale with exp(shift))
          kinetic_data.(fn{it}).mean(ind_mean)    = my_mean * exp(my_shift(ind_mean));
          kinetic_data.(fn{it}).std(ind_mean)     = my_std  * exp(my_shift(ind_mean));
          kinetic_data.(fn{it}).median(ind_mean)  = exp(my_mean_ln + my_shift(ind_mean));
          kinetic_data.(fn{it}).mean_ln(ind_mean) = my_mean_ln + my_shift(ind_mean);
          kinetic_data.(fn{it}).std_ln(ind_mean)  = my_std_ln;
      end

    end

  end

end

if ~isempty(network)

  if isfield(kinetic_data,'KM')
    kinetic_data.KM.lower(network.N'==0) = nan;
    kinetic_data.KM.upper(network.N'==0) = nan;
  end

  if isfield(kinetic_data,'KA')
    kinetic_data.KA.lower(network.regulation_matrix <=0) = nan;
    kinetic_data.KA.upper(network.regulation_matrix <=0) = nan;
  end

  if isfield(kinetic_data,'KI')
    kinetic_data.KI.lower(network.regulation_matrix >=0) = nan;
    kinetic_data.KI.upper(network.regulation_matrix >=0) = nan;
  end

end


% -----------------------------------------------------------------

function q = complete_log_bounds(q)

% Where a bound is given only on non-logarithmic scale (lower, upper),
% insert the corresponding value on logarithmic scale (lower_ln, upper_ln)

for bound = {'lower','upper'}
  b    = bound{1};
  b_ln = [b '_ln'];
  if isfield(q,b)
    if ~isfield(q,b_ln) || isempty(q.(b_ln)),
      q.(b_ln) = nan * q.(b);
    end
    ind = find(isfinite(q.(b)) & ~isfinite(q.(b_ln)) & q.(b) > 0);
    q.(b_ln)(ind) = log(q.(b)(ind));
  end
end
