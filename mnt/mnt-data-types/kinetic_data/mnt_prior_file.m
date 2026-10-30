function prior_file = mnt_prior_file()

% prior_file = mnt_prior_file()
%
% Return the filename of the default prior table (SBtab format)
%
% Other packages can use their own prior tables, by passing the filename to
% 'parameter_balancing_prior'

prior_file = fullfile(fileparts(mfilename('fullpath')), 'mnt_prior.tsv');
