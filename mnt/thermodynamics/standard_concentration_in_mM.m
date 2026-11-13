function [c, label] = standard_concentration_in_mM(text)

% [c, label] = standard_concentration_in_mM(text)
%
% Standard concentration, given as text, in mM
%
% text:  '1 M', '1M', 'M' (result 1000) or '1 mM', '1mM', 'mM' (result 1); other values give an error
% label: '1M' or '1mM' (form used for the table attribute StandardConcentration in files written by MNT)
%
% Standard concentrations are declared in SBtab tables by the table attribute StandardConcentration
% (for equilibrium constants, standard chemical potentials and standard reaction Gibbs free energies).
% Matlab functions in MNT use 1 mM internally. See also standard_concentration_convert

switch strrep(text, ' ', '')
  case {'M', '1M'},   c = 1000; label = '1M';
  case {'mM', '1mM'}, c = 1;    label = '1mM';
  otherwise, error(sprintf('Unknown standard concentration ''%s'' (allowed: ''1 M'' or ''1 mM'')', text));
end
