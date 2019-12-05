function pb_dir = pb_DIR()
  
pb_dir = [fileparts(which(mfilename)) filesep '..' filesep];
