function dir = get_pb_directory()
  
dir = [fileparts(which(mfilename)) filesep '..'];
