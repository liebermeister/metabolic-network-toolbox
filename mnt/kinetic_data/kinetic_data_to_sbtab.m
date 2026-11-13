function [sbtab_table, M] = kinetic_data_to_sbtab(kinetic_data, network, table_id, standard_concentration)

% KINETIC_DATA_TO_SBTAB Convert a kinetic data struct (see kinetic_data_load) into an SBtab table object
%
% [sbtab_table, M] = kinetic_data_to_sbtab(kinetic_data, network, table_id, standard_concentration)
%
% kinetic_data:  kinetic data struct (e.g. from kinetic_data_load)
% network:       network struct
% table_id:      (optional) TableID of the table, default 'ParameterData'
% standard_concentration: (optional) standard concentration ('1 mM' or '1 M') for equilibrium constants,
%                standard chemical potentials and standard reaction Gibbs free energies in the table; default '1 mM'
%
% sbtab_table:   SBtab table object (TableType 'Quantity'; columns !QuantityType, !Compound, !Reaction,
%                !Unit, !Median, !Mean, !Std, !GeomStd, !Lower, !Upper); it can be read by kinetic_data_load
% M:             the same table as a cell array of strings (including the declaration line)
%
% In the kinetic_data struct, equilibrium constants and standard Gibbs free energies refer to the
% standard concentration 1 mM. In the table, they are converted to the standard concentration given
% (water is not counted), and the table declares it with the attribute StandardConcentration.
%
% To save the table to a file, use kinetic_data_save. See also: kinetic_data_load

eval(default('table_id','''ParameterData''','standard_concentration','''1 mM'''));
if isempty(table_id), table_id = 'ParameterData'; end
if isempty(standard_concentration), standard_concentration = '1 mM'; end

[~, standard_concentration_label] = standard_concentration_in_mM(standard_concentration);
kinetic_data = kinetic_data_convert_standard_concentration(kinetic_data, network, standard_concentration);

T = kinetic_data_print(kinetic_data,network,0);

M = [];

attribute_line_written = 0;

for it = 1:length(T)
  if size(T{it},1)>1
    if attribute_line_written == 0;
      M = [M; T{it}(1:end,:)];
      attribute_line_written = 1;
    else
      M = [M; T{it}(2:end,:)];
    end
  end
end

% add attribute line

M = [[{sprintf('!!SBtab TableID=''%s'' TableType=''Quantity'' StandardConcentration=''%s''', table_id, standard_concentration_label)}, repmat({''},1,size(M,2)-1)]; M];

sbtab_table = sbtab_table_load([], M);


% ---------------------------------------------------------------

function kinetic_data = kinetic_data_convert_standard_concentration(kinetic_data, network, standard_concentration)
% convert Keq, dmu0 and mu0 from 1 mM to the standard concentration given

if standard_concentration_in_mM(standard_concentration) == 1, return; end

conversions = {'Keq', 'ln Keq', length(network.actions); 'dmu0', 'dmu0', length(network.actions); 'mu0', 'mu0', length(network.metabolites)};

for it = 1:size(conversions,1)
  symbol = conversions{it,1};
  if ~isfield(kinetic_data, symbol), continue; end
  shift = standard_concentration_convert(zeros(conversions{it,3},1), conversions{it,2}, network, '1 mM', standard_concentration);
  x = kinetic_data.(symbol);
  if strcmp(symbol, 'Keq')
    % multiplicative quantity: values, means, standard deviations and bounds scale with exp(shift);
    % geometric standard deviations and std of ln values do not change
    for fn = {'median', 'mean', 'std', 'lower', 'upper'}
      if isfield(x, fn{1}), x.(fn{1}) = x.(fn{1}) .* exp(shift); end
    end
    for fn = {'mean_ln', 'lower_ln', 'upper_ln'}
      if isfield(x, fn{1}), x.(fn{1}) = x.(fn{1}) + shift; end
    end
  else
    % additive quantity: values, means and bounds are shifted; standard deviations do not change
    for fn = {'median', 'mean', 'lower', 'upper'}
      if isfield(x, fn{1}), x.(fn{1}) = x.(fn{1}) + shift; end
    end
  end
  kinetic_data.(symbol) = x;
end
