function bigg_match_metabolite_ids_display(report, title_text)

% BIGG_MATCH_METABOLITE_IDS_DISPLAY Print the report of bigg_match_metabolite_ids
%
% bigg_match_metabolite_ids_display(report, title_text)
%
% report      report struct from bigg_match_metabolite_ids
% title_text  (optional) first line, e.g. 'Mapping of metabolite data to the model'

if ~exist('title_text', 'var') || isempty(title_text), title_text = 'Mapping of metabolite ids (data -> model)'; end

fprintf('%s: %d of %d data metabolites mapped\n', title_text, report.n_mapped, report.n_data);
show('  not found in the model (ignored)', report.not_found);
show('  ambiguous: several model metabolites (ignored)', report.ambiguous);
show('  duplicate: model metabolite already mapped (ignored)', report.duplicate);
if ~isempty(report.compartment_ignored)
  pairs = strcat(report.compartment_ignored, {' -> '}, report.compartment_ignored_model_ids);
  show('  mapped with different compartment', pairs);
end
if ~isempty(report.external_data_ids)
  pairs = strcat(report.external_data_ids, {' -> '}, report.external_model_ids);
  show('  mapped to external metabolites of the model', pairs);
end


function show(label, ids)

if ~isempty(ids)
  fprintf('%s: %s\n', label, strjoin(ids(:)', ', '));
end
