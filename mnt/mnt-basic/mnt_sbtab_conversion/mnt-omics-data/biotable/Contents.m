% Matlab data structure for omics data matrices (metabolites, fluxes, proteins, expression)
%
% REPRESENTATION IN MATLAB: CELL STRUCTURE
%
% The first n fields contain names, IDs etc. of the biological elements (as string lists)
% The number n of these fields may vary
%
% data_table.SampleName: list of sample names (strings)
% data_table.SampleTime: vector of sample times (numbers)
% data_table.DataMean:   matrix of data mean values
% data_table.DataStd:    matrix of data standard deviations
% data_table.Info:       N x 2 cellstruct with strings (additional information)
%
%
% REPRESENTATION IN FILES: TAB-SEPARATED FILES (THREE POSSIBLE FORMATS)
%
% (I) Simple tab-separated files (two files are necessary)
%
%   1. File [name].tsv: contains everything except for .Info
%    First n columns: contents of the first n fields
%    Remaining columns: data, columns appear in the order
%        "[sample1] Mean" "[sample1] Std" "[sample2] Mean" "[sample2] Std" etc
%    First row: column headers
%    Second row: sample times
%    Remaining rows: contents
%
%   2. File [name]_info.tsv: contains the contents of .Info
%
% (II) SBtab files
%
% (III) OpenBIS files
%
%
% FUNCTIONS
%
% biotable_construct_timeseries     - construct a biotable
% biotable_data_construct           - construct biotables for metabolites, fluxes, proteins and transcripts
% biotable_read, biotable_write     - read and write biotable files
% biotable_read_multiple, biotable_write_multiple - same, for metabolites, fluxes, proteins and transcripts
% biotable_choose_items, biotable_remove, biotable_remove_nan_items, biotable_choose_samples - select data
% biotable_align                    - reorder the elements
% biotable_concatenate, biotable_combine_experiments - concatenate samples
% biotable_join_replicates          - collect data from several replicates
% biotable_combine_replicates       - combine replicates by replicate regression (requires replicate-regression)
% biotable_replicate_mean_std       - mean values and standard deviations of replicates
% biotable_interpolate, biotable_time_interpolation - interpolation in time
% biotable_network_interpolate      - fill in missing values from neighbours in the network
% biotable_smooth_mean, biotable_omit_outliers - preprocessing
% biotable_add_initial_time_point, biotable_invent_standard_deviations, biotable_clean - other changes
% biotable_data_consistency         - compare flux data with metabolite and expression data
% openBIS_BE_default_pars           - default options for writing files in OpenBIS format
