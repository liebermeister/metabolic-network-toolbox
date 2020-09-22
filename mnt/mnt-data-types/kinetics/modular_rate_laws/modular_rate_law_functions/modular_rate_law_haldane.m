function modular_rate_law_haldane(network)
  
[~,~,Keq] = modular_kcat_to_KV_Keq(network.N,network.kinetics,network.kinetics.Kcatf,network.kinetics.Kcatr);

if norm(log(network.kinetics.Keq ./ Keq)) > 10^-10,
  display('Equilibrium constants (from kinetics | computed from other parameters):');
  [network.kinetics.Keq,Keq]
  error('Haldane relationship is violated!')
else
  display('The Haldane relationships seem to be satisfied')
end
