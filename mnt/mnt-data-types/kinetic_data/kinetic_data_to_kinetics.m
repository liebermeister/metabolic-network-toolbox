function kinetics = kinetic_data_to_kinetics(kinetics, kinetic_data, network);

% [kinetics, kinetics_std_nat] = kinetic_data_to_kinetics(kinetics, kinetic_data);
%
% Generate 'kinetics' data structure from median values in 'kinetic_data' data structure
%
% the field .std_nat contains the standard deviations of quantities in the "natural scaling" (linear or logarithmic, depending on the type of quantity
  

if isempty(kinetics), 
  kinetics = set_kinetics(network,'cs');
end

fn = union(fieldnames(kinetics),fieldnames(kinetic_data));

for it = 1:length(fn),
  if isfield(kinetic_data,fn{it}),
    kinetics.(fn{it}) = kinetic_data.(fn{it}).median;
  end
end

if isfield(kinetics,'dmu0'),
  if find(isnan(kinetics.Keq) .* isfinite(kinetics.dmu0)),
    display('kinetic_data_to_kinetics: inserting Keq values from known standard Delta G values');
  end
  kinetics.Keq(isnan(kinetics.Keq)) = exp(-1/RT * kinetics.dmu0(isnan(kinetics.Keq)));
end

if ~isfield(kinetics,'Kcatf'), 
  if isfield(kinetic_data,'Kcatf'),
    kinetics.Kcatf = kinetic_data.Kcatf.median;
    kinetics.Kcatr = kinetic_data.Kcatr.median;
  else
  [kinetics.Kcatf, kinetics.Kcatr] = modular_KV_Keq_to_kcat(network.N,kinetics,kinetics.KV,kinetics.Keq,kinetics.KM,kinetics.h);
  end
end

% ------------------------
% kinetics_std

kinetics_std_nat = struct;

for it = 1:length(fn),
  if isfield(kinetic_data,fn{it}),
    if isfield(kinetic_data.(fn{it}),'std_ln'),
      kinetics_std_nat.(fn{it}) = kinetic_data.(fn{it}).std_ln;
    else
      kinetics_std_nat.(fn{it}) = kinetic_data.(fn{it}).std;
    end
  end
end

if isfield(kinetic_data,'Kcatf'),
  kinetics_std_nat.Kcatf = kinetic_data.Kcatf.std_ln;
  kinetics_std_nat.Kcatr = kinetic_data.Kcatr.std_ln;
end

kinetics.STD_nat = kinetics_std_nat;