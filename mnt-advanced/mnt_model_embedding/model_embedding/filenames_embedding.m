function filenames = filenames_embedding(model_name)

filenames.network_dir     = ['/home/wolfram/matlab/projects/model_embedding/models/' model_name '/'];
filenames.network_file    = [filenames.network_dir model_name];
filenames.table_positions = [filenames.network_dir model_name '_Position.tsv'];
filenames.psfile_dir      = ['/home/wolfram/projekte/dfg-projekt_dynamics_and_function/projekte/model_embedding/ps-files/model_embedding/' model_name];
filenames.result_file = [filenames.network_dir '/' model_name '_result'];

