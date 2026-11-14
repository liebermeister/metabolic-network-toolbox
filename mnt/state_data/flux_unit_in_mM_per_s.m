function f = flux_unit_in_mM_per_s(unit, context)

% f = flux_unit_in_mM_per_s(unit, context)
%
% Conversion factor from a flux unit to mM/s (value in mM/s = f * value in 'unit')
%
% unit:    text '<concentration unit>/<time unit>', e.g. 'mM/s', 'uM/min', 'mmol/l/h'
%          (concentration units as in concentration_unit_in_mM; time units s, min, h)
% context: optional text for error messages (e.g. file name and quantity)
%
% Unknown units give an error. Matlab functions in MNT use mM/s for fluxes.
% See also concentration_unit_in_mM

eval(default('context', ''''''));

u = strrep(char(unit), ' ', '');
k = find(u == '/', 1, 'last');
if isempty(k)
  error_unit(unit, context);
end
switch u(k+1:end)
  case {'s', 'sec', 'second'},  t = 1;
  case {'min', 'minute'},       t = 60;
  case {'h', 'hour'},           t = 3600;
  otherwise, error_unit(unit, context);
end
f = concentration_unit_in_mM(u(1:k-1), context) / t;


% -------------------------------------------------------------------

function error_unit(unit, context)

if ~isempty(context), context = [' (' context ')']; end
error('Unknown flux unit ''%s''%s; allowed: <concentration unit>/<time unit>, e.g. mM/s, uM/min, mM/h', char(unit), context);
