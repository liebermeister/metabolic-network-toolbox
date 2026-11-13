function [sbtab_table, M] = kinetic_data_to_sbtab(kinetic_data, network, table_id)

% KINETIC_DATA_TO_SBTAB Convert a kinetic data struct (see kinetic_data_load) into an SBtab table object
%
% [sbtab_table, M] = kinetic_data_to_sbtab(kinetic_data, network, table_id)
%
% kinetic_data:  kinetic data struct (e.g. from kinetic_data_load)
% network:       network struct
% table_id:      (optional) TableID of the table, default 'ParameterData'
%
% sbtab_table:   SBtab table object (TableType 'Quantity'; columns !QuantityType, !Compound, !Reaction,
%                !Unit, !Median, !Mean, !Std, !GeomStd, !Lower, !Upper); it can be read by kinetic_data_load
% M:             the same table as a cell array of strings (including the declaration line)
%
% Equilibrium constants and standard Gibbs free energies refer to the standard concentration 1 mM
% (as in kinetic_data); the table declares this with the attribute StandardConcentration='1mM'.
%
% To save the table to a file, use kinetic_data_save. See also: kinetic_data_load

eval(default('table_id','''ParameterData'''));

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

M = [[{sprintf('!!SBtab TableID=''%s'' TableType=''Quantity'' StandardConcentration=''1mM''', table_id)}, repmat({''},1,size(M,2)-1)]; M];

sbtab_table = sbtab_table_load([], M);
