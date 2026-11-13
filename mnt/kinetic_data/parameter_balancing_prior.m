function [parameter_prior, parameter_prior_filename] = parameter_balancing_prior(omit_quantities, parameter_prior_filename, verbose)

% parameter_prior = parameter_balancing_prior(omit_quantities,parameter_prior_filename)
%
% Data structure containing information about parameter priors (read from the Prior table file)
%
% parameter_prior_filename: prior table (file name, SBtab table object, or SBtab document object
%                           containing a table of type 'QuantityInfo'); by default, the table given by 'mnt_prior_file'
%
% The standard concentration of the prior table (table attribute StandardConcentration, '1 mM' or '1 M')
% is stored in parameter_prior.standard_concentration. It refers to the prior values of equilibrium
% constants, standard chemical potentials and standard reaction Gibbs free energies. A prior table that
% does not declare it is read as '1 mM' (with a warning). Functions that use the prior values for a
% network (pb_task, kinetic_data_complete, model balancing) convert them to 1 mM, the standard concentration
% used internally in MNT (see prior_standard_concentration_shift).

if ~exist('sbtab_version','file'),
  error('For this function, the SBtab toolbox must be installed');
end

eval(default('omit_quantities','[]','parameter_prior_filename','[]','verbose','0'));

if isempty(parameter_prior_filename),
  parameter_prior_filename = mnt_prior_file;
end

% The prior file contains a single table; its TableID differs between versions ('Prior', 'PriorTable')
if ischar(parameter_prior_filename),
  parameter_prior_sbtab = sbtab_table_load(parameter_prior_filename);
else
  [parameter_prior_sbtab, parameter_prior_filename] = sbtab_table_from_input(parameter_prior_filename, {'PriorTable','Prior'}, 'QuantityInfo');
end

if verbose,
  display(sprintf('o Using parameter prior file %s', parameter_prior_filename));
end
parameter_prior       = sbtab_table_get_all_columns(parameter_prior_sbtab);

if isfield(parameter_prior,'PriorGeometricStd'),
  %display('Using geometric standard deviations');
  my_std = log(cell_string2num(parameter_prior.PriorGeometricStd)) / log(10);
  ind = isfinite(my_std);
  parameter_prior.PriorStd(ind) = num2cellstr(my_std(ind));
end

% -----------------------------------------------------------------------------
% omit certain quantities from the prior table

keep = 1:length(parameter_prior.QuantityType);

if length(omit_quantities),
  keep = setdiff(keep,label_names(omit_quantities,parameter_prior.QuantityType));
end

fn = fieldnames(parameter_prior);
for it = 1:length(fn),
  parameter_prior.(fn{it}) = parameter_prior.(fn{it})(keep);
end

for it = 1:length(parameter_prior.Symbol),
  parameter_prior.symbol_index.(parameter_prior.Symbol{it})=it;
end

% standard concentration of thermodynamic quantities in the prior table
table_attributes = sbtab_table_get_attributes(parameter_prior_sbtab);
if isfield(table_attributes,'StandardConcentration')
  parameter_prior.standard_concentration = table_attributes.StandardConcentration;
else
  parameter_prior.standard_concentration = '1 mM';
  warning('The prior table does not declare a standard concentration (table attribute StandardConcentration); I assume 1 mM');
end
standard_concentration_in_mM(parameter_prior.standard_concentration); % check the value
