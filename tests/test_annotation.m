function tests = test_annotation

% Tests for the BiGG-KEGG id conversion (folder mnt/annotation)
%
% Run:  results = runtests('test_annotation')

tests = functiontests(localfunctions);


function test_bigg_to_kegg(testCase)

[m, r] = convert_bigg_to_kegg({'glc[e]'; 'atp'; 'pyr[c]'; 'nonexistent'}, {'PGI'; 'PYK'; 'nonexistent'});
verifyEqual(testCase, m, {'C00031'; 'C00002'; 'C00022'; ''});
verifyEqual(testCase, r(2), {'R00200'});
verifyEqual(testCase, r(3), {''});


function test_kegg_to_bigg(testCase)

[m, r] = convert_kegg_to_bigg({'C00002'; 'C00022'; 'C99999'}, {'R00200'; 'R99999'});
verifyEqual(testCase, m, {'atp'; 'pyr'; ''});
verifyEqual(testCase, r, {'PYK'; ''});


function test_extra_mappings(testCase)

filename = [tempname '.tsv'];
fid = fopen(filename, 'w');
fprintf(fid, '!!!SBtab Document=''test''\n\n!!SBtab TableID=''Compound'' TableType=''Compound''\n!Identifiers:kegg.compound\t!Identifiers:bigg.metabolite\nC00002\tmy_atp\nC99999\tnew\n');
fclose(fid);
testCase.addTeardown(@() delete(filename));
m = convert_kegg_to_bigg({'C00002'; 'C99999'; 'C00022'}, {}, struct('extra_mapping_files', {{filename}}));
verifyEqual(testCase, m, {'my_atp'; 'new'; 'pyr'});
m = convert_kegg_to_bigg({'C00002'; 'C00022'}, {}, struct('extra_mapping_files', {{filename}}, 'metabolites_from_main_table', 0));
verifyEqual(testCase, m, {'my_atp'; ''});
