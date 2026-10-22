function setup_test_path()

% setup_test_path()
%
% Add the toolbox to the Matlab path and check that matlab-utils is installed
% (called at the start of each test file)

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'mnt')));

if ~exist('lp_solve','file'),
  error('Please add matlab-utils (https://github.com/liebermeister/matlab-utils) to the Matlab path');
end
