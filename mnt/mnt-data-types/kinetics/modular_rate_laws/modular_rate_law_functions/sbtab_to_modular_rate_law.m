function [kinetics, sbtab_table, other_parameters] = sbtab_to_modular_rate_law(network, file_kinetic_data, options)

% [kinetics, sbtab_table, other_parameters] = sbtab_to_modular_rate_law(network,file_kinetic_data,options)
%
% Build 'kinetics' data structure from kinetic data (in SBtab file)
%
% Please note:
%   Keq values are assumed to refer to a standard concentration of mM.
%   If in the data file, Keq values refer to a standard concentration of M, 
%   this needs to be indicated by the table attribute StandardConcentration='M' (in the Parameter table)
%   The values will then be converted by this function.
%
% Similar functions: 'ms_import_kinetic', 'sbtab_to_modular_rate_law_via_kinetic_data'

eval(default('options','struct'));

options_default = struct('use_sbml_ids',1,'kinetic_law','cs','verbose','0');
options         = join_struct(options_default,options);

switch options.kinetic_law,
  case {'cs','ms','rp','ma','fm'}, % UPDATE rate law names!
     kinetics = set_kinetics(network,options.kinetic_law);
  otherwise, error('Conversion is only possible for modular rate laws');
end

if isstr(file_kinetic_data),
  sbtab_table = sbtab_table_load(file_kinetic_data);
else
  %% assume that file_kinetic_data contains already an sbtab data structure
  sbtab_table = file_kinetic_data;
end

table_attributes = sbtab_table_get_attributes(sbtab_table);
if isfield(table_attributes,'StandardConcentration'),
  standard_concentration = table_attributes.StandardConcentration;
else
  standard_concentration = 'mM';
end

QuantityType = sbtab_table_get_column(sbtab_table,'QuantityType');
if sbtab_table_has_column(sbtab_table,'Value')
  Value        = cell_string2num(sbtab_table_get_column(sbtab_table,'Value'));
else
  Value        = cell_string2num(sbtab_table_get_column(sbtab_table,'Mode'));
end
Compound     = sbtab_table_get_column(sbtab_table,'Compound');
Reaction     = sbtab_table_get_column(sbtab_table,'Reaction');
compound_ind = label_names(Compound,network.metabolites);
reaction_ind = label_names(Reaction,network.actions);
sdm = setdiff(network.metabolites,Compound);
sdr = setdiff(network.actions,Reaction);
if length(sdm)+length(sdr),
  sdm
  sdr
  error('Information about metabolites and / or reactions missing in data file');
end
[nr,nm,nx,ind_KM,ind_KA,ind_KI,nKM,nKA,nKI] = network_numbers(network);

ind_Keq = find(strcmp(QuantityType,'equilibrium constant'));
ind_Keq = ind_Keq(label_names(Reaction(ind_Keq),network.actions));
kinetics.Keq = Value(ind_Keq);

if length(ind_Keq),
  switch standard_concentration,
    case {'M','1M','1 M'},
      %% convert equilibrium constants: standard concentration M -> standard concentration mM
      ind_water = network_find_water(network);
      N_for_dissolved_compounds = network.N;
      N_for_dissolved_compounds(ind_water,:) = 0;
      molecularity_mismatch = full(sum(N_for_dissolved_compounds,1));
      kinetics.Keq = kinetics.Keq .* 1000.^column(molecularity_mismatch);
      display('Converting equilibrium constants from a convention with standard concentration M -> standard concentration mM.');
    case {'mM','1mM','1 mM'},
    otherwise
      warning('I assume that equilibrium constants given in the data file refer to a standard concentration of mM. If this is not the case, please indicate this by adding the table attribute StandardConcentration=''M'' to your data file');
  end
end

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

ind_CompoundMass = find(strcmp(QuantityType,'molecular mass'));
ll = label_names(network.metabolites,Compound(ind_CompoundMass));
other_parameters.metabolite_mass     = nan * ones(nm,1);
other_parameters.metabolite_mass(find(ll)) = Value(ind_CompoundMass(ll(find(ll))));

ind_EnzymeMass = find(strcmp(QuantityType,'protein molecular mass'));
ll = label_names(network.actions,Reaction(ind_EnzymeMass));
other_parameters.enzyme_mass   = nan * ones(nr,1);
other_parameters.enzyme_mass(find(ll))     = Value(ind_EnzymeMass(ll(find(ll))));

if isempty(kinetics.Keq),
  error('No Keq values found');
end

% -----------------------
% Double check consistency between kinetic constants

modular_rate_law_haldane(network,kinetics);
modular_rate_law_haldane_kcat(network,kinetics);
