% ---------------------------------------------------------------------------------
% Smaller version of E. coli model (Noor et al 2016) containing only glycolysis
% ---------------------------------------------------------------------------------

network = network_sbml_import('ecoli_noor_2016.xml','/home/wolfram/matlab/wolf_packages/mnt/mnt/mnt-extensions/mnt_parameter_balancing/parameter_balancing/../models/',0);

glycolysis_reactions = {'PTS_RPTSsy','PGI_R02740','PFK_R04779','ALD_R01070','TIM_R01015','GAP_R01061','PGK_R01512','PGM_R01518','PGH_R00658','PYK_R00200','PDH_R00209'};
 
ind_r = label_names(glycolysis_reactions, network.actions);
subnetwork = network_subnetwork(network,[],ind_r,1);

network_sbml_export(subnetwork, 0, 'E. coli Noor (2016), glycolysis subnetwork', '/home/wolfram/matlab/wolf_packages/mnt/mnt/mnt-extensions/mnt_parameter_balancing/models/ecoli_noor_2016_glycolysis/ecoli_noor_2016_glycolysis.xml');
  
network_to_sbtab(subnetwork, struct('filename','/home/wolfram/matlab/wolf_packages/mnt/mnt/mnt-extensions/mnt_parameter_balancing/models/ecoli_noor_2016_glycolysis/ecoli_noor_2016_glycolysis.tsv'));
  
