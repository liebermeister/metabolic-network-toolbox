function check_flux_stationarity(network,v)

% check_flux_stationarity(network,v)

ind_int = find(network.external==0);
mismatch = norm(network.N(ind_int,:) * v)
if mismatch,
  network.metabolites(ind_int(find(mismatch)))
end
stationary = [mismatch<10^-10]