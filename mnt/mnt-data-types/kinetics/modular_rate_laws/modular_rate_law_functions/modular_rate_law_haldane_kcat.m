function modular_rate_law_haldane_kcat(network,kinetics)
  
eval(default('kinetics','network.kinetics'));
  
% Double check consistency between kinetic constants

[computed_Kcatf, computed_Kcatr] = modular_KV_Keq_to_kcat(network.N,kinetics,kinetics.KV,kinetics.Keq,kinetics.KM,kinetics.h);

ind_violf = find(abs(log(computed_Kcatf) - log(kinetics.Kcatf)) > 0.01);
ind_violr = find(abs(log(computed_Kcatr) - log(kinetics.Kcatr)) > 0.01);
if length(ind_violf) + length(ind_violr), 
  display('- modular_rate_law_haldane_kcat.m: Haldane relationships seem to be violated (mismatch in forward Kcat values)');
  reactions = network.actions(ind_violf);
  Kcatf_existent_and_required = pm([kinetics.Kcatf(ind_violf),computed_Kcatf(ind_violf)],reactions)
  reactions = network.actions(ind_violr);
  kcatr_in_model = kinetics.Kcatr(ind_violr);
  Kcatr_existent_and_required = pm([kinetics.Kcatr(ind_violr),computed_Kcatr(ind_violr)],reactions)
end 
