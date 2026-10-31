function sbtab = sbtab_network_state_save(network,c,v,filename,dmu,options, keq)

% sbtab = sbtab_network_state_save(network,c,v,filename,dmu,options, keq)
%
% Save a 'network' data structure, concentration vector 'c', and flux vector 'v' to SBtab file

eval(default('options','struct'));
options = join_struct(struct('verbose',1), options);

eval(default('dmu','[]','keq','[]'));

sbtab = network_to_sbtab(network,struct('c',c,'v',v,'dmu',dmu,'keq',keq));

sbtab_document_save(sbtab,filename,0,options.verbose);
