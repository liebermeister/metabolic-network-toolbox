function d = mnt_BASEDIR()

% d = mnt_BASEDIR()
%
% Main folder of the Metabolic Network Toolbox

d = [fileparts(which(mfilename)) filesep '..' filesep];
