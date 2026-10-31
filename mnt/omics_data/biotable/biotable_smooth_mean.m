function data = biotable_smooth_mean(data)

% data = biotable_smooth_mean(data)
%
% preprocess data table in order to reach a smooth behavior
% of the geometric mean as a function of time

mm            = median(log(data.DataMean),1,'omitnan');
mm_smooth     = my_bfilt(my_bfilt(mm));
data.DataMean = exp( log(data.DataMean) + repmat(mm_smooth-mm,size(data.DataMean,1),1));
