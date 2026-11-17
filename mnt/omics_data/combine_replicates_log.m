function [ln_avg, ln_std] = combine_replicates_log(ln_rep_avg, ln_rep_std, ln_std_default)

% COMBINE_REPLICATES_LOG Combine replicate measurements on log scale (mean and error bars)
%
% [ln_avg, ln_std] = combine_replicates_log(ln_rep_avg, ln_rep_std, ln_std_default)
%
% Combine replicates into a single data set with error bars, on log scale (assuming gaussian errors).
%
% ln_rep_avg:     log values, array (e.g. metabolites x time points x replicates; replicates in dimension 3);
%                 missing values: NaN
% ln_rep_std:     standard deviations of these log values (same size; empty if the data have no error bars)
% ln_std_default: standard deviation for data points without error bars
%
% With error bars:    weighted mean with weights 1/ln_rep_std.^2, and the standard deviation of this mean
%                     (missing error bars are replaced by ln_std_default)
% Without error bars: mean, and the standard error of the mean (standard deviation between replicates / sqrt(n));
%                     if this cannot be computed, for instance for one replicate, ln_std_default is used
%
% (Until 2026-11-18 the IMF function imf_combine_replicates.)
%
% See also: biotable_replicate_mean_std

is_measured = isfinite(ln_rep_avg);

if ~isempty(ln_rep_std)
  ln_rep_std(is_measured & ~(ln_rep_std>0)) = ln_std_default;
  weights = is_measured ./ ln_rep_std.^2;
  weights(~is_measured) = 0;
  ln_rep_avg(~is_measured) = 0;
  ln_avg = sum(weights .* ln_rep_avg, 3) ./ sum(weights, 3);
  ln_std = 1 ./ sqrt(sum(weights, 3));
  ln_std(sum(weights, 3)==0) = nan;
else
  n_measured = sum(is_measured, 3);
  ln_avg = mean(ln_rep_avg, 3, 'omitnan');
  ln_std = std(ln_rep_avg, 0, 3, 'omitnan') ./ sqrt(n_measured);
  ln_std(isfinite(ln_avg) & ~(ln_std>0)) = ln_std_default;
end
