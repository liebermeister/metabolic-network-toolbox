% Tests for the biotable data structure (omics data)

setup_test_path;
data = biotable_construct_timeseries('MetaboliteName', {'A';'B'}, 0:5, [1 2 3 4 5 6; 1 nan 3 nan 5 6]);

%% construct
assert(isequal(data.MetaboliteName, {'A';'B'}));
assert(isequal(data.SampleTime, (0:5)'));
assert(all(isnan(data.DataStd(:))));

%% choose and remove items and samples
assert(isequal(biotable_choose_items(data, 2).MetaboliteName, {'B'}));
assert(isequal(biotable_remove(data, 1).MetaboliteName, {'B'}));
assert(isequal(biotable_choose_samples(data, [1 3]).SampleTime, [0;2]));

%% interpolation in time
data2 = biotable_time_interpolation(data, [0.5; 1.5]);
assert(isequal(data2.DataMean(1,:), [1.5 2.5]));
data3 = biotable_interpolate(data, [0.5 4.5]);
assert(isequal(data3.DataMean(1,:), [1.5 5.5]));

%% smoothing
data4 = biotable_smooth_mean(biotable_construct_timeseries('MetaboliteName', {'A'}, 0:5, exp(1:6)));
assert(isequal(size(data4.DataMean), [1 6]));

%% omit outliers
data5 = biotable_construct_timeseries('MetaboliteName', {'A'}, 0:5, [1 1 1 100 1 1]);
data5 = biotable_omit_outliers(data5, 0.2, 0);
assert(isnan(data5.DataMean(4)));
