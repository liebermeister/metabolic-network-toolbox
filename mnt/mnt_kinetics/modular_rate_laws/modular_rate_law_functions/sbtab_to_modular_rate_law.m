function [kinetics, sbtab_table] = sbtab_to_modular_rate_law(network, file_kinetic_data, options)

% [kinetics, sbtab_table] = sbtab_to_modular_rate_law(network,file_kinetic_data,options)

eval(default('options','struct'));

options_default = struct('use_sbml_ids',1,'kinetic_law','cs','verbose','0');
options         = join_struct(options_default,options);

switch options.kinetic_law,
  case {'cs','ms','rp','ma','fm'}, % UPDATE rate law names!
     kinetics = set_kinetics(network,options.kinetic_law);
  otherwise, error('Conversion is only possible for modular rate law');
end

sbtab_table  = sbtab_table_load(file_kinetic_data);
QuantityType = sbtab_table_get_column(sbtab_table,'QuantityType');
Value        = cell_string2num(sbtab_table_get_column(sbtab_table,'Value'));
Compound     = sbtab_table_get_column(sbtab_table,'Compound');
Reaction     = sbtab_table_get_column(sbtab_table,'Reaction');
compound_ind = label_names(Compound,network.metabolites);
reaction_ind = label_names(Reaction,network.actions);

[nr,nm,nx,ind_KM,ind_KA,ind_KI,nKM,nKA,nKI] = network_numbers(network);

ind_Keq = find(strcmp(QuantityType,'equilibrium constant'));
ind_Keq = ind_Keq(label_names(Reaction(ind_Keq),network.actions));
kinetics.Keq = Value(ind_Keq);

ind_KV = find(strcmp(QuantityType,'catalytic rate constant geometric mean'));
ind_KV = ind_KV(label_names(Reaction(ind_KV),network.actions));
kinetics.KV = Value(ind_KV);

ind_c = find(strcmp(QuantityType,'concentration'));
ind_c = ind_c(label_names(Compound(ind_c),network.metabolites));
kinetics.c = Value(ind_c);

ind_u = find(strcmp(QuantityType,'concentration of enzyme'));
ind_u = ind_u(label_names(Reaction(ind_u),network.actions));
kinetics.u = Value(ind_u);

ind_KM  = find(strcmp(QuantityType,'Michaelis constant'));
ind_KMc = label_names(Compound(ind_KM),network.metabolites);
ind_KMr = label_names(Reaction(ind_KM),network.actions);
kinetics.KM = zeros(nr,nm);
kinetics.KM(sub2ind([nr,nm],ind_KMr,ind_KMc)) = Value(ind_KM);

ind_KI  = find(strcmp(QuantityType,'inhibitory constant'));
ind_KIc = label_names(Compound(ind_KI),network.metabolites);
ind_KIr = label_names(Reaction(ind_KI),network.actions);
kinetics.KI = zeros(nr,nm);
kinetics.KI(sub2ind([nr,nm],ind_KIr,ind_KIc)) = Value(ind_KI);

ind_KA  = find(strcmp(QuantityType,'activation constant'));
ind_KAc = label_names(Compound(ind_KA),network.metabolites);
ind_KAr = label_names(Reaction(ind_KA),network.actions);
kinetics.KA = zeros(nr,nm);
kinetics.KA(sub2ind([nr,nm],ind_KAr,ind_KAc)) = Value(ind_KA);

[kinetics.Kcatf, kinetics.Kcatr] = modular_KV_Keq_to_kcat(network.N,kinetics,kinetics.KV,kinetics.Keq,kinetics.KM,kinetics.h);
