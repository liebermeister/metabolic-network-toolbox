function [data_table, report] = sbtab_orient_data_table(network, data_table, options)

% [data_table, report] = sbtab_orient_data_table(network, data_table, options)
%
% Convert the values in an SBtab data table to the reaction orientation of a model.
%
% For each table row, the orientation of the reaction (same as in the model, or opposite) is
% determined by sbtab_reaction_orientation (from reaction formulas given in the table, or from
% the table attribute ReactionOrientation). In rows with the opposite orientation, the values are
% converted according to their quantity type (table QuantityOrientationDependence, file
% quantity_orientation_dependence.tsv in the folder of this function):
%   sign:     values change their sign (fluxes, Gibbs free energies of reaction); Min and Max are exchanged
%   inverse:  values are inverted (equilibrium constants); Min and Max are exchanged
%   swap:     the quantity type is replaced (e.g. substrate catalytic rate constant -> product
%             catalytic rate constant)
%   none:     no change
% Then, the reaction formulas in the converted rows are reversed, and the table attribute
% ReactionOrientation is set to options.model_orientation (if all rows referring to model reactions
% have a known orientation).
%
% Inputs
%   network      network structure (see sbtab_reaction_orientation)
%   data_table   SBtab table object
%   options      struct with (optional) fields
%     .match_data_by      see sbtab_reaction_orientation (default 'ModelElementId')
%     .model_orientation  see sbtab_reaction_orientation
%     .dependence_file    file with table QuantityOrientationDependence
%                         (default: quantity_orientation_dependence.tsv in the folder of this function)
%     .verbose            (default 0)
%
% Outputs
%   data_table   converted table
%   report       struct with fields
%     .reversed_reactions      reaction IDs (table entries) of the converted rows
%     .reversed_quantities     quantity types of the converted rows
%     .unknown_orientation     reaction IDs of rows with orientation-dependent quantities whose
%                              orientation is unknown (values are not converted)
%     .unknown_quantities      quantity types without a known orientation dependence (in reversed rows)

eval(default('options','struct'));
options_default = struct('match_data_by', 'ModelElementId', 'model_orientation', [], ...
                         'dependence_file', [], 'verbose', 0);
options = join_struct(options_default, options);

if isempty(options.dependence_file)
  options.dependence_file = fullfile(fileparts(mfilename('fullpath')), 'quantity_orientation_dependence.tsv');
end

report = struct('reversed_reactions', {{}}, 'reversed_quantities', {{}}, ...
                'unknown_orientation', {{}}, 'unknown_quantities', {{}});

n_rows = sbtab_table_n_rows(data_table);
if n_rows == 0, return; end

[orientation, ind_reaction] = sbtab_reaction_orientation(network, data_table, options);

% -------------------------------------
% quantity types of the rows

quantity_types = sbtab_table_get_column(data_table, 'QuantityType', 0, 0);
if isempty(quantity_types)
  if isfield(data_table.attributes, 'QuantityType')
    quantity_types = repmat({data_table.attributes.QuantityType}, n_rows, 1);
  else
    quantity_types = repmat({''}, n_rows, 1);
  end
end
quantity_types = column(quantity_types);

dependence_table = sbtab_table_load(options.dependence_file);
dep_quantity = sbtab_table_get_column(dependence_table, 'QuantityType');
dep_type     = sbtab_table_get_column(dependence_table, 'ReactionOrientationDependence');
dep_swap     = sbtab_table_get_column(dependence_table, 'SwapWith');

row_dependence = repmat({'unknown'}, n_rows, 1);
row_swap       = repmat({''}, n_rows, 1);
for it = 1:n_rows
  ll = find(strcmp(quantity_types{it}, dep_quantity));
  if length(ll)
    row_dependence{it} = dep_type{ll(1)};
    row_swap{it}       = dep_swap{ll(1)};
  end
end

reaction_ids = sbtab_table_get_column(data_table, 'Reaction', 0, 0);
if isempty(reaction_ids), reaction_ids = repmat({''}, n_rows, 1); end

% -------------------------------------
% columns to be converted

[names, kinds] = value_columns(data_table);

% -------------------------------------
% convert rows with opposite orientation

ind_reverse = find(orientation == -1);

for it = column(ind_reverse)'
  switch row_dependence{it}
    case 'none',
      continue
    case 'sign',
      data_table = convert_row(data_table, it, names, kinds, @(x) -x);
    case 'inverse',
      data_table = convert_row(data_table, it, names, kinds, @(x) 1./x);
    case 'swap',
      quantity_types{it} = row_swap{it};
    otherwise,
      report.unknown_quantities = [report.unknown_quantities; quantity_types(it)];
      continue
  end
  report.reversed_reactions  = [report.reversed_reactions;  reaction_ids(it)];
  report.reversed_quantities = [report.reversed_quantities; quantity_types(it)];
end

if sbtab_table_has_column(data_table, 'QuantityType')
  data_table = sbtab_table_set_column(data_table, 'QuantityType', quantity_types);
end

% reverse the formulas in the converted rows (all rows with opposite orientation)

for formula_column = {'ReactionFormula', 'ReactionFormula:Identifiers:kegg.compound'}
  f = sbtab_table_get_column(data_table, formula_column{1}, 0, 0);
  if ~isempty(f)
    for it = column(ind_reverse)'
      f{it} = reverse_formula(f{it});
    end
    data_table = sbtab_table_set_column(data_table, formula_column{1}, f);
  end
end

% rows with unknown orientation (orientation-dependent quantities only)

is_dependent  = ~strcmp(row_dependence, 'none');
ind_unknown   = find(isfinite(ind_reaction) & isnan(orientation) & is_dependent);
report.unknown_orientation = reaction_ids(ind_unknown);

if ~isempty(options.model_orientation) && isempty(ind_unknown)
  data_table.attributes.ReactionOrientation = options.model_orientation;
end

if options.verbose
  table_name = '';
  if isfield(data_table.attributes, 'TableID'), table_name = data_table.attributes.TableID; end
  if length(report.reversed_reactions)
    display(sprintf('  Table %s: values converted to the model orientation in %d rows (reactions %s)', ...
                    table_name, length(report.reversed_reactions), strjoin(unique(report.reversed_reactions)', ', ')));
  end
  if length(ind_unknown)
    warning(sprintf('Table %s: unknown reaction orientation in %d rows (reactions %s); values not converted', ...
                    table_name, length(ind_unknown), strjoin(unique(report.unknown_orientation)', ', ')));
  end
  if length(report.unknown_quantities)
    warning(sprintf('Table %s: orientation dependence unknown for quantities %s; values not converted', ...
                    table_name, strjoin(unique(report.unknown_quantities)', ', ')));
  end
end


% ------------------------------------------------------------------

function [names, kinds] = value_columns(data_table)

% numerical columns that depend on the orientation
% kinds: 'value' (converted), 'min', 'max' (converted and exchanged), 'std' (see convert_row)

names = {}; kinds = {};
all_names = [column(data_table.column.column_names); column(data_table.uncontrolled.headers)];
all_names = all_names(cellfun(@ischar, all_names));
all_names = all_names(~cellfun(@isempty, all_names));
for it = 1:length(all_names)
  nm = regexprep(all_names{it}, '^!', '');
  if sum(strcmp(nm, {'Value', 'Mean', 'Median', 'Mode'})) || ~isempty(regexp(nm, '_(Mean|Median|Value)$', 'once'))
    names{end+1} = nm; kinds{end+1} = 'value';
  elseif strcmp(nm, 'Min') || ~isempty(regexp(nm, '_Min$', 'once'))
    names{end+1} = nm; kinds{end+1} = 'min';
  elseif strcmp(nm, 'Max') || ~isempty(regexp(nm, '_Max$', 'once'))
    names{end+1} = nm; kinds{end+1} = 'max';
  elseif strcmp(nm, 'Std') || ~isempty(regexp(nm, '_Std$', 'once'))
    names{end+1} = nm; kinds{end+1} = 'std';
  elseif isempty(regexp(nm, '[:]', 'once')) && is_state_column(data_table, nm)
    names{end+1} = nm; kinds{end+1} = 'value';  % e.g. columns of a QuantityMatrix (one per state)
  end
end


% ------------------------------------------------------------------

function flag = is_state_column(data_table, nm)

% columns without '!' (QuantityMatrix states), except annotation columns
flag = sum(strcmp(nm, data_table.uncontrolled.headers)) && ...
       ~sum(strcmpi(nm, {'Comment', 'Reference', 'Organism', 'Unit', 'Name', 'direction'}));


% ------------------------------------------------------------------

function data_table = convert_row(data_table, it, names, kinds, fun)

% apply 'fun' to the values of row 'it'; Min and Max are converted and exchanged;
% Std columns: for x -> -x, the Std does not change; for x -> 1/x, Std -> Std / x^2 (first order)

values = cell(size(names));
for k = 1:length(names)
  c = sbtab_table_get_column(data_table, names{k}, 0, 0);
  values{k} = c;
end

is_inverse = isequal(fun(2), 0.5);
x_value = nan;
ind_value = find(strcmp(kinds, 'value'));
if length(ind_value)
  x_value = str2double(values{ind_value(1)}{it});
end

ind_min = find(strcmp(kinds, 'min'));
ind_max = find(strcmp(kinds, 'max'));

for k = 1:length(names)
  x = str2double(values{k}{it});
  if isnan(x), continue; end
  switch kinds{k}
    case 'value',
      values{k}{it} = num2str(fun(x), 8);
    case 'std',
      if is_inverse && isfinite(x_value) && x_value ~= 0
        values{k}{it} = num2str(x / x_value^2, 8);
      end
  end
end

% Min and Max: convert and exchange (pairs with the same prefix)
for k = column(ind_min)'
  partner = find(strcmp(names, regexprep(names{k}, 'Min$', 'Max')));
  if isempty(partner), continue; end
  x_min = str2double(values{k}{it}); x_max = str2double(values{partner}{it});
  new_min = fun(x_max); new_max = fun(x_min);
  values{k}{it}       = num2str_nan(new_min);
  values{partner}{it} = num2str_nan(new_max);
end

for k = 1:length(names)
  data_table = sbtab_table_set_column(data_table, names{k}, values{k});
end


% ------------------------------------------------------------------

function s = num2str_nan(x)

if isnan(x), s = 'NaN'; else, s = num2str(x, 8); end


% ------------------------------------------------------------------

function f = reverse_formula(f)

if isempty(f), return; end
[tok, sep] = regexp(f, '\s*(<=>|<->|->|=)\s*', 'split', 'match');
if length(tok) == 2
  f = [strtrim(tok{2}) ' <=> ' strtrim(tok{1})];
end
