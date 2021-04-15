function modular_rate_law_haldane(network,kinetics,fignum)

eval(default('kinetics','[]','fignum','9999'));

if isempty(kinetics), kinetics = network.kinetics; end

network.kinetics = kinetics;
  
% modular_rate_law_haldane(network)

[~,~,Keq_from_Haldane] = modular_kcat_to_KV_Keq(network.N, network.kinetics);

epsilon = 0.001;

if max(abs(log10(network.kinetics.Keq ./ Keq_from_Haldane))) > epsilon,
  [dum,ind] = sort(-abs(log10([network.kinetics.Keq./Keq_from_Haldane])));
  ind = ind(dum<-epsilon);
  figure(fignum);
  plot(log10(network.kinetics.Keq),log10(Keq_from_Haldane),'.');
  text(log10(network.kinetics.Keq(ind)), log10(Keq_from_Haldane(ind)), strrep(network.actions(ind),'_','-'));
  xlabel('log10 Keq in model'); ylabel('log10 Keq obtained from Haldane relationship');
  title('- (modular_rate_law_haldane.m): Some Haldane relationships are violated!')
  display(sprintf('- modular_rate_law_haldane.m: Some Haldane relationships are violated! (difference of log10 Keq values > %f)',epsilon))
else
  %display('All Haldane relationships are satisfied')
end
