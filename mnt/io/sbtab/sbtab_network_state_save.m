function sbtab = sbtab_network_state_save(network,c,v,filename,dmu,options, keq)

% sbtab = sbtab_network_state_save(network,c,v,filename,dmu,options, keq)
%
% Save a 'network' data structure, concentration vector 'c', and flux vector 'v' to SBtab file
% (optionally, also reaction Gibbs free energies 'dmu' and equilibrium constants 'keq')
%
% keq refers to the standard concentration 1 mM. Option standard_concentration ('1 mM' (default) or
% '1 M'): standard concentration used in the file (declared by the table attribute StandardConcentration)

eval(default('options','struct'));
options = join_struct(struct('verbose',1,'standard_concentration','1 mM'), options);

eval(default('dmu','[]','keq','[]'));

sbtab = network_to_sbtab(network,struct('c',c,'v',v,'dmu',dmu,'keq',keq,'standard_concentration',options.standard_concentration));

sbtab_document_save(sbtab,filename,0,options.verbose);
