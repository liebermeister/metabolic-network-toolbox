function [kinetics,S_ext] = par2standard(network,par,n_ext)

% [kinetics,S_ext] = par2standard(network,par,n_ext)

if ~exist('n_ext','var'), n_ext = 0; end

kinetics = network.kinetics;
kinetics.type = 'standard';
z=1;
for it=1:length(kinetics.reactions),
 for itt = 1:length(kinetics.reactions{it}.parameters),
   ll = kinetics.reactions{it}.sizes(itt);
   kinetics.reactions{it}=setfield(kinetics.reactions{it},kinetics.reactions{it}.parameters{itt},par(z:z+ll-1));
   z=z+ll;
 end
end

S_ext = par(end-n_ext+1:end);