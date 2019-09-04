model_file   = [pb_DIR '/models/ecoli_noor_2016.xml'];
data_file    = [pb_DIR '/models/ecoli_noor_2016_data.tsv'];
prior_file   = [pb_DIR '/config/pb_prior.tsv'];
options_file = [pb_DIR '/config/pb_options.tsv'];
%output_file  = '/tmp/ecoli_noor_2016_result';
output_file  = [];
model_name   = 'E. coli central metabolism, Noor et al (2016)';
result       = parameter_balancing(model_file, output_file, data_file, prior_file, options_file, model_name,struct('flag_check',1,'preferred_data_element_ids','id'));
