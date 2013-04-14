function isStationary = check_flux_stationarity(network,v)

% check_flux_stationarity(network,v)

ind_int  = find(network.external==0);
mismatch = network.N(ind_int,:) * v;
mismatch = mismatch(find(abs(mismatch)>10^-10));
if norm(mismatch),
  display('Internal metabolites not balanced:');
  network.metabolites(ind_int(find(mismatch)))
  mismatch
end
isStationary = [norm(mismatch)<10^-10];