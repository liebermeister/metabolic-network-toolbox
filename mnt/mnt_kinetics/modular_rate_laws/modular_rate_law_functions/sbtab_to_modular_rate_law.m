function [kinetics, sbtab_table, other_parameters] = sbtab_to_modular_rate_law(network, file_kinetic_data, options)

% [kinetics, sbtab_table, other_parameters] = sbtab_to_modular_rate_law(network,file_kinetic_data,options)

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

ind_Kcatf = find(strcmp(QuantityType,'substrate catalytic rate constant'));
ind_Kcatf = ind_Kcatf(label_names(Reaction(ind_Kcatf),network.actions));
kinetics.Kcatf = Value(ind_Kcatf);

ind_Kcatr = find(strcmp(QuantityType,'product catalytic rate constant'));
ind_Kcatr = ind_Kcatr(label_names(Reaction(ind_Kcatr),network.actions));
kinetics.Kcatr = Value(ind_Kcatr);

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

ind_c = find(strcmp(QuantityType,'concentration'));
ll    = label_names(Compound(ind_c),network.metabolites);
ind_c = ind_c(find(ll));
kinetics.c = Value(ind_c);

ind_u = find(strcmp(QuantityType,'concentration of enzyme'));
ll    = label_names(Reaction(ind_u),network.actions);
ind_u = ind_u(find(ll));
kinetics.u = Value(ind_u);

ind_CompoundMass = find(strcmp(QuantityType,'compound molecular weight'));
ll = label_names(network.metabolites,Compound(ind_CompoundMass));
other_parameters.compound_mass     = nan * ones(nm,1);
other_parameters.compound_mass(find(ll)) = Value(ind_CompoundMass(ll(find(ll))));

ind_EnzymeMass = find(strcmp(QuantityType,'enzyme molecular weight'));
ll = label_names(network.actions,Reaction(ind_EnzymeMass));
other_parameters.enzyme_mass   = nan * ones(nr,1);
other_parameters.enzyme_mass(find(ll))     = Value(ind_EnzymeMass(ll(find(ll))));

[computed_Kcatf, computed_Kcatr] = modular_KV_Keq_to_kcat(network.N,kinetics,kinetics.KV,kinetics.Keq,kinetics.KM,kinetics.h);

if norm(log(computed_Kcatf) - log(kinetics.Kcatf)) > 0.01 * length(computed_Kcatf), 
  warning('Given Kcat values and Kcat values computed from other parameters do not match'); 
end 

if norm(log(computed_Kcatr) - log(kinetics.Kcatr)) > 0.01 * length(computed_Kcatr), 
  warning('Given Kcat values and Kcat values computed from other parameters do not match'); 
end 
