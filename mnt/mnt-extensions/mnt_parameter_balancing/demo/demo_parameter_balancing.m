% --------------------------------------------------
% Parameter balancing demo script
% Usage example for matlab function 'parameter_balancing'
%
% Example model: E. coli central carbon metabolism, from Noor et al. 2016 (PLoS Comp Biol)
%
% Input files
%   model_file:   (SBML format)  from models/ecoli_noor_2016.xml      
%   data_file :   (SBtab format) from models/ecoli_noor_2016_data.tsv 
%   prior_file:   default from subdirectory config/
%   options_file: default from subdirectory config/
%
% Results are written to tmp directory
% --------------------------------------------------

model_name   = 'E. coli central metabolism, Noor et al (2016)';
model_file   = [pb_DIR '/models/ecoli_noor_2016.xml'];
data_file    = [pb_DIR '/models/ecoli_noor_2016_data.tsv'];
prior_file   = [pb_DIR '/config/pb_prior.tsv'];
options_file = [pb_DIR '/config/pb_options.tsv'];
output_file  = [tempdir 'demo_pb_ecoli_noor_2016.tsv'];

[result, pb_options] = parameter_balancing(model_file, output_file, data_file, prior_file, options_file, model_name,struct('flag_check',1,'preferred_data_element_ids','id'));
