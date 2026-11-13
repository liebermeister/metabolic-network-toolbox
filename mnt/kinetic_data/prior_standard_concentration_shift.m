function shift = prior_standard_concentration_shift(parameter_prior, symbol, network)

% shift = prior_standard_concentration_shift(parameter_prior, symbol, network)
%
% Shift that converts prior values (median, mean, bounds) of a quantity from the standard
% concentration of the prior table to 1 mM, the standard concentration used internally in MNT
%
% parameter_prior  prior structure (see parameter_balancing_prior), with field .standard_concentration
%                  (if the field is missing, 1 mM is assumed)
% symbol           symbol of the quantity in the prior table (e.g. 'Keq', 'dmu0', 'mu0', 'KM')
% network          network structure (needed only if the prior table does not use 1 mM)
%
% shift            values to be added to the prior values in their natural scaling:
%                  'Keq': vector (one value per reaction) for ln Keq (multiply Keq by exp(shift))
%                  'dmu0': vector (one value per reaction) in kJ/mol
%                  'mu0': vector (one value per metabolite) in kJ/mol (0 for water)
%                  other quantities: 0 (they do not depend on the standard concentration)

shift = 0;

if ~isfield(parameter_prior, 'standard_concentration')
  return
end

if standard_concentration_in_mM(parameter_prior.standard_concentration) == 1
  return
end

switch symbol
  case {'Keq', 'dmu0', 'mu0'}
    if isempty(network)
      error('The prior table uses the standard concentration %s; to convert the prior values for %s to 1 mM, a network is needed', parameter_prior.standard_concentration, symbol);
    end
end

switch symbol
  case 'Keq'
    shift = standard_concentration_convert(zeros(length(network.actions),1), 'ln Keq', network, parameter_prior.standard_concentration, '1 mM');
  case 'dmu0'
    shift = standard_concentration_convert(zeros(length(network.actions),1), 'dmu0', network, parameter_prior.standard_concentration, '1 mM');
  case 'mu0'
    shift = standard_concentration_convert(zeros(length(network.metabolites),1), 'mu0', network, parameter_prior.standard_concentration, '1 mM');
end
