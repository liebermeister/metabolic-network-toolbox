function f = sbtab_table_concentration_factor(sbtab_table, table_label)

% f = sbtab_table_concentration_factor(sbtab_table, table_label)
%
% Conversion factor from the concentration unit of an SBtab table to mM (values in mM = f * values in the table)
%
% The unit is taken from the table attribute Unit or, if missing, from the column Unit (all rows must have the
% same unit). If neither is given, mM is assumed (f = 1). Unknown units and different units in one table give
% an error. table_label (optional): name used in error messages.
%
% See also concentration_unit_in_mM

eval(default('table_label', '''SBtab table'''));

f = 1;
if isempty(sbtab_table), return; end

unit = '';
attributes = sbtab_table.attributes;
if isfield(attributes, 'Unit') && ~isempty(attributes.Unit)
  unit = attributes.Unit;
elseif sbtab_table_has_column(sbtab_table, 'Unit')
  units = sbtab_table_get_column(sbtab_table, 'Unit');
  units = unique(units(~cellfun(@isempty, strtrim(units))));
  if length(units) > 1
    error('Different units in %s: %s', table_label, strjoin(column(units)', ', '));
  end
  if length(units), unit = units{1}; end
end

if ~isempty(unit)
  f = concentration_unit_in_mM(unit, table_label);
end
