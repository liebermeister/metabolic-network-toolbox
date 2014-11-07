function network = network_build_from_sum_formulae_list(reactions)

% network = network_build_from_sum_formulae_list(reactions)
%
% wrapper around network_build_from_sum_formulae for matlab reaction list as input

dum = struct; 
dum.SumFormula = reactions;
network = network_build_from_sum_formulae([],[],dum);
