function d = mnt_basedir()

% MNT_BASEDIR Return the main directory of the Metabolic Network Toolbox repository
%
% d = mnt_basedir()

d = fileparts(fileparts(fileparts(which(mfilename))));
