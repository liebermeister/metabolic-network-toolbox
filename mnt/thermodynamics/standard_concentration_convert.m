function values = standard_concentration_convert(values, quantity, network, from, to, temperature)

% values = standard_concentration_convert(values, quantity, network, from, to, temperature)
%
% Convert thermodynamic quantities from one standard concentration to another
%
% values      vector (one value per reaction, or per metabolite for 'mu0') or matrix with one
%             row per reaction or metabolite (e.g. several metabolic states); NaN values stay NaN
% quantity    'Keq'    equilibrium constants (also for lower and upper bounds, medians and
%                      geometric means; not for arithmetic standard deviations)
%             'ln Keq' natural logarithms of equilibrium constants
%             'dmu0'   standard reaction Gibbs free energies (kJ/mol)
%             'mu0'    standard chemical potentials (kJ/mol)
%             'dmu'    reaction Gibbs free energies (kJ/mol): do not depend on the standard
%                      concentration and are returned unchanged
% network     network structure (fields 'metabolites', 'N' and, optionally, annotations used
%             to recognise water)
% from, to    standard concentrations as text ('1 M', '1 mM', ...; see standard_concentration_in_mM)
% temperature temperature in K (default 298.15), used for 'dmu0' and 'mu0'
%
% With c0 = standard concentration and m = sum of stoichiometric coefficients of a reaction:
%   Keq(to)  = Keq(from) * (c0_from/c0_to)^m
%   dmu0(to) = dmu0(from) - RT * m * ln(c0_from/c0_to)
%   mu0(to)  = mu0(from)  - RT * ln(c0_from/c0_to)
% Water (recognised by network_find_water) is not counted, since its activity is set to 1.
% Matlab functions in MNT use 1 mM internally.

eval(default('temperature','298.15'));

L = log(standard_concentration_in_mM(from) / standard_concentration_in_mM(to));

if L == 0 || strcmp(quantity, 'dmu') || isempty(values)
  return
end

ind_water = network_find_water(network);

switch quantity
  case {'Keq', 'ln Keq', 'dmu0'}
    N = network.N;
    N(ind_water,:) = 0;
    m = full(column(sum(N,1)));
    if size(values,1) ~= length(m)
      error('Number of values (%d) differs from number of reactions (%d)', size(values,1), length(m));
    end
    switch quantity
      case 'Keq',    values = values .* repmat(exp(m * L), 1, size(values,2));
      case 'ln Keq', values = values  + repmat(m * L, 1, size(values,2));
      case 'dmu0',   values = values  - repmat(RT(temperature) * m * L, 1, size(values,2));
    end
  case 'mu0'
    if size(values,1) ~= length(network.metabolites)
      error('Number of values (%d) differs from number of metabolites (%d)', size(values,1), length(network.metabolites));
    end
    shift = - RT(temperature) * L * ones(length(network.metabolites),1);
    shift(ind_water) = 0;
    values = values + repmat(shift, 1, size(values,2));
  otherwise
    error('Unknown quantity ''%s''', quantity);
end
