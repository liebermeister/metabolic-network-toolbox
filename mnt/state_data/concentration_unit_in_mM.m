function f = concentration_unit_in_mM(unit, context)

% f = concentration_unit_in_mM(unit, context)
%
% Conversion factor from a concentration unit to mM (value in mM = f * value in 'unit')
%
% unit:    text, e.g. 'M', 'mM', 'uM', 'µM', 'nM', 'mol/l', 'mmol/l', 'mol/L', 'mmol/L', 'umol/l', 'nmol/l'
%          (spaces are ignored)
% context: optional text for error messages (e.g. file name and quantity)
%
% Unknown units give an error. Matlab functions in MNT use mM for all concentrations
% (metabolite and enzyme concentrations, KM, KA, KI values). See also standard_concentration_in_mM

eval(default('context', ''''''));

u = strrep(strrep(char(unit), ' ', ''), 'µ', 'u');
u = strrep(u, 'L', 'l');

switch u
  case {'M',  'mol/l'},   f = 1000;
  case {'mM', 'mmol/l'},  f = 1;
  case {'uM', 'umol/l'},  f = 0.001;
  case {'nM', 'nmol/l'},  f = 10^-6;
  otherwise
    if ~isempty(context), context = [' (' context ')']; end
    error('Unknown concentration unit ''%s''%s; allowed: M, mM, uM, nM (or mol/l, mmol/l, umol/l, nmol/l)', char(unit), context);
end
