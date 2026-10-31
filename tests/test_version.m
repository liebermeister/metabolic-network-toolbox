% Tests for the version number

setup_test_path;

%% The version in mnt_version and in Contents.m agree
contents_file = fullfile(fileparts(mfilename('fullpath')), '..', 'mnt', 'Contents.m');
assert(contains(fileread(contents_file), ['% Version ' mnt_version]));
