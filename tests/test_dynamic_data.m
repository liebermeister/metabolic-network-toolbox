% Tests for the conversion between dynamic data (omics data) and SBtab

setup_test_path;

%% omics data: dynamic data -> SBtab document -> file -> dynamic data
data = struct('MetaboliteName', {{'A';'B'}}, 'SampleName', {{'t1';'t2';'t3'}}, 'SampleTime', [0;1;2], ...
              'DataMean', [1 2 3; 4 5 6], 'DataStd', [0.1 0.2 0.3; 0.4 0.5 0.6], 'Info', {{'Organism','E. coli'}});
d = dynamic_data_2_sbtab(data, 'Test');
filename = [tempname '.tsv'];
sbtab_document_save(d, filename, 0, 0);
data2 = sbtab_2_dynamic_data(filename);
delete(filename);
assert(isequal(data2.MetaboliteName, data.MetaboliteName));
assert(isequal(data2.SampleName(:), data.SampleName));
assert(isequal(data2.SampleTime(:), data.SampleTime));
assert(isequal(data2.DataMean, data.DataMean));
assert(isequal(data2.DataStd, data.DataStd));
assert(isequal(data2.Info, data.Info));
