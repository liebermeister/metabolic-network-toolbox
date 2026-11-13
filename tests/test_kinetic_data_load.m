% Tests for kinetic_data_load:
% - a table with a Mean column: rows that give only a Median are kept
% - Keq data are clipped to the prior bounds after the conversion of the standard concentration

setup_test_path;
if ~exist('sbtab_version','file'), return; end

% R1: A + B -> C (sum of stoichiometric coefficients m = -1); R2: C -> D (m = 0)
N = [-1 0; -1 0; 1 -1; 0 1];
network = network_construct(N, [1;1], [], {'A';'B';'C';'D'}, {'R1';'R2'});

parameter_prior = parameter_balancing_prior;
parameter_prior.standard_concentration = '1M';
ind_keq = find(strcmp(parameter_prior.Symbol,'Keq'));
keq_max = eval(parameter_prior.UpperBound{ind_keq});
keq_min = eval(parameter_prior.LowerBound{ind_keq});
unit_keq = parameter_prior.Unit{ind_keq};
unit_km  = parameter_prior.Unit{strcmp(parameter_prior.Symbol,'KM')};

% data in 1 mM: R1 Keq = keq_min / 100 is inside the bounds (in 1 mM, the lower bound of R1 is
% keq_min / 1000); R2 Keq = 10 * keq_max is outside the bounds (m = 0, no conversion)
lines = {"!!SBtab TableID='ParameterData' TableType='Quantity' StandardConcentration='1mM'", ...
         sprintf('!QuantityType\t!Reaction\t!Compound\t!Unit\t!Median\t!Mean\t!Std'), ...
         sprintf('equilibrium constant\tR1\t\t%s\t%g\t\t', unit_keq, keq_min / 100), ...
         sprintf('equilibrium constant\tR2\t\t%s\t%g\t\t', unit_keq, keq_max * 10), ...
         sprintf('Michaelis constant\tR1\tA\t%s\t\t0.5\t0.1', unit_km), ...
         sprintf('Michaelis constant\tR1\tB\t%s\t0.2\t\t', unit_km)};
filename = [tempname '.tsv'];
fid = fopen(filename, 'w'); fprintf(fid, '%s\n', lines{:}); fclose(fid);

opt = struct('verbose', 0, 'use_sbml_ids', 0, 'use_kegg_ids', 0, 'enforce_ranges', 1);
[~, kinetic_data] = evalc('kinetic_data_load({''equilibrium constant'',''Michaelis constant''}, parameter_prior, network, filename, opt)');
delete(filename);

%% Mean and Median in one table
assert(abs(kinetic_data.KM.median(1,2) - 0.2) < 1e-12);   % row with only a Median
assert(abs(kinetic_data.KM.mean(1,1) - 0.5) < 1e-12);     % row with a Mean
assert(isfinite(kinetic_data.KM.mean_ln(1,2)) && isfinite(kinetic_data.KM.mean_ln(1,1)));
assert(abs(kinetic_data.KM.mean_ln(1,2) - log(0.2)) < 1e-12);

%% Keq: converted first (here 1 mM -> 1 mM, no change), then clipped with the bounds converted to 1 mM
assert(abs(kinetic_data.Keq.median(1) / (keq_min / 100) - 1) < 1e-10);   % not clipped
assert(abs(kinetic_data.Keq.median(2) / keq_max - 1) < 1e-10);           % clipped
assert(abs(kinetic_data.Keq.mean_ln(2) - log(keq_max)) < 1e-10);

%% a table with Keq data but without StandardConcentration is rejected; without Keq data it is accepted
lines_keq = {"!!SBtab TableID='ParameterData' TableType='Quantity'", ...
             sprintf('!QuantityType\t!Reaction\t!Compound\t!Unit\t!Median'), ...
             sprintf('equilibrium constant\tR1\t\t%s\t2', unit_keq)};
lines_km  = {"!!SBtab TableID='ParameterData' TableType='Quantity'", ...
             sprintf('!QuantityType\t!Reaction\t!Compound\t!Unit\t!Median'), ...
             sprintf('Michaelis constant\tR1\tA\t%s\t0.3', unit_km)};
filename = [tempname '.tsv'];
fid = fopen(filename, 'w'); fprintf(fid, '%s\n', lines_keq{:}); fclose(fid);
failed = 0;
try, evalc('kinetic_data_load({''equilibrium constant''}, parameter_prior, network, filename, opt)'); catch, failed = 1; end
assert(failed);
fid = fopen(filename, 'w'); fprintf(fid, '%s\n', lines_km{:}); fclose(fid);
[~, kinetic_data] = evalc('kinetic_data_load({''equilibrium constant'',''Michaelis constant''}, parameter_prior, network, filename, opt)');
delete(filename);
assert(abs(kinetic_data.KM.median(1,1) - 0.3) < 1e-12);
