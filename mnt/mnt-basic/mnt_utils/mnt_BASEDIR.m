function d = mnt_BASEDIR()

d = [fileparts(which(mfilename)) filesep '..' filesep '..' filesep ];
