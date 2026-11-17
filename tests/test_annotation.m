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


function test_strip_metabolite_ids(testCase)

% only suffixes are removed (until 2026-11-15, substrings were removed anywhere)
s = bigg_kegg_strip_metabolite_ids({'glc__D_e'; '13dpg_h'; 'glc[e]'; '13_cis_retn_c'; 'mal__L'; 'acon_C'; 'my_atp'});
verifyEqual(testCase, s, {'glc'; '13dpg'; 'glc'; '13_cis_retn'; 'mal'; 'acon_C'; 'my_atp'});
[b, c] = bigg_metabolite_id_base({'icit[c]'; 'glc__D_e'; 'mal__L'; '13dpg_h'; 'pyr'});
verifyEqual(testCase, b, {'icit'; 'glc__D'; 'mal__L'; '13dpg'; 'pyr'});
verifyEqual(testCase, c, {'c'; 'e'; ''; 'h'; ''});
verifyEqual(testCase, convert_bigg_to_kegg({'13dpg_h'; '13dpg[c]'}), {'C00236'; 'C00236'});


function test_match_metabolite_ids(testCase)

model_ids = {'icit[c]'; 'glc[c]'; 'glc[e]'; '13dpg[c]'; 'mal__L[c]'; 'pi[e]'};
external  = [0; 0; 1; 0; 0; 1];
data_ids  = {'icit'; 'glc'; 'glc_e'; '13dpg_h'; 'mal'; 'pi'; 'xyz'; 'icit_c'};
[ind, report] = bigg_match_metabolite_ids(model_ids, data_ids, external);
verifyEqual(testCase, ind, [1; 2; 3; 4; 0; 6]);
verifyEqual(testCase, report.not_found, {'mal'; 'xyz'});
verifyEqual(testCase, report.compartment_ignored, {'13dpg_h'; 'pi'});   % pi -> pi[e]: no internal pi
verifyEqual(testCase, report.duplicate, {'icit_c'});
verifyEqual(testCase, report.external_data_ids, {'glc_e'; 'pi'});
verifyTrue(testCase, report.has_problems);

%% convert_kegg_to_bigg: empty KEGG ids give empty BiGG ids
[m, r] = convert_kegg_to_bigg({'C00022'; ''}, {''; 'R00200'});
assert(isempty(m{2}) && isempty(r{1}));
assert(strcmp(m{1}, 'pyr') && ~isempty(r{2}));
