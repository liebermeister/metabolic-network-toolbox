function [theta_along_flux,A,A_along_flux] = driving_force(c,v,network,keq)
  
% theta     = -Delta G/RT; in flux direction
% A_along_flux = -Delta G; in flux direction
% A         = -Delta G; direction by convention

if ~exist('keq','var'), keq = network.kinetics.Keq; end
  
dmu = - RT * [log(keq) - network.N' * log(c)];
A = - dmu;
A_along_flux = - sign(v) .* dmu;
theta_along_flux     = A_along_flux / RT;

if sum(theta_along_flux(v~=0)<0), 
  [theta_along_flux, v]
  error('Thermodynamic flux constraint violated'); 
end 
