function [ind, report] = bigg_match_metabolite_ids(model_ids, data_ids, model_external)

% BIGG_MATCH_METABOLITE_IDS Match metabolite ids of a data set to the metabolite ids of a model
%
% [ind, report] = bigg_match_metabolite_ids(model_ids, data_ids, model_external)
%
% model_ids       cell array of model metabolite ids (e.g. BiGG ids 'icit', 'icit[c]', 'glc[e]', 'icit_c')
% data_ids        cell array of metabolite ids in the data set
% model_external  (optional) boolean vector: external metabolites of the model (default: none)
%
% ind     vector (length of model_ids): index of the matching data row for each model metabolite (0: none)
% report  struct describing the mapping (see below; shown by bigg_match_metabolite_ids_display)
%
% Compartments are written either as labels in brackets ('glc[e]') or as BiGG suffixes ('glc_e',
% '13dpg_h'); see bigg_metabolite_id_base. Each data id is assigned to (at most) one model metabolite:
%   1. identical id
%   2. same base id and compatible compartment ('' and 'c' (cytosol) count as the same compartment),
%      e.g. 'icit' -> 'icit[c]', 'glc_e' -> 'glc[e]'
%   3. same base id, compartment ignored (e.g. data '13dpg_h' -> model '13dpg[c]'); reported
% If several model metabolites qualify in a step, internal metabolites are preferred; if there is still
% more than one candidate, the data id is not mapped (reported as ambiguous). If several data ids map to
% the same model metabolite, the first one is used (reported). Stereo suffixes ('__L', '__D') are part of
% the id and are not ignored (data 'mal__L' does not match model 'mal').
%
% report fields
%   .n_data, .n_mapped                 numbers of data ids and of mapped data ids
%   .mapped_data_ids, .mapped_model_ids        all mapped pairs
%   .not_found                         data ids without a model metabolite
%   .ambiguous                         data ids with several equally good model metabolites
%   .compartment_ignored               data ids mapped in step 3 (with .compartment_ignored_model_ids)
%   .duplicate                         data ids not used because another data id maps to the same metabolite
%   .external_data_ids, .external_model_ids    data ids mapped to external model metabolites
%   .has_problems                      true if any of not_found, ambiguous, compartment_ignored, duplicate
%                                      is non-empty

model_ids = model_ids(:);
data_ids  = data_ids(:);
if ~exist('model_external', 'var') || isempty(model_external)
  model_external = false(size(model_ids));
end
model_external = logical(model_external(:));

[model_base, model_comp] = bigg_metabolite_id_base(model_ids);
[data_base,  data_comp]  = bigg_metabolite_id_base(data_ids);
model_comp(strcmp(model_comp, 'c')) = {''};
data_comp(strcmp(data_comp, 'c'))   = {''};

ind = zeros(length(model_ids), 1);
report = struct('n_data', length(data_ids), 'n_mapped', 0, 'mapped_data_ids', {{}}, 'mapped_model_ids', {{}}, ...
                'not_found', {{}}, 'ambiguous', {{}}, 'compartment_ignored', {{}}, 'compartment_ignored_model_ids', {{}}, ...
                'duplicate', {{}}, 'external_data_ids', {{}}, 'external_model_ids', {{}}, 'has_problems', false);

for j = 1:length(data_ids)
  if ~ischar(data_ids{j}) || isempty(data_ids{j}), continue; end
  same_base = strcmp(model_base, data_base{j});
  steps = {strcmp(model_ids, data_ids{j}), same_base & strcmp(model_comp, data_comp{j}), same_base};
  i = []; step = 0;
  for it = 1:3
    candidates = find(steps{it});
    if length(candidates) > 1 && any(~model_external(candidates))
      candidates = candidates(~model_external(candidates));
    end
    if length(candidates) == 1
      i = candidates; step = it; break
    elseif length(candidates) > 1
      report.ambiguous{end+1, 1} = data_ids{j}; step = -1; break
    end
  end
  if step == -1, continue; end
  if isempty(i)
    report.not_found{end+1, 1} = data_ids{j};
    continue
  end
  if ind(i) ~= 0
    report.duplicate{end+1, 1} = data_ids{j};
    continue
  end
  ind(i) = j;
  report.mapped_data_ids{end+1, 1}  = data_ids{j};
  report.mapped_model_ids{end+1, 1} = model_ids{i};
  if step == 3
    report.compartment_ignored{end+1, 1}           = data_ids{j};
    report.compartment_ignored_model_ids{end+1, 1} = model_ids{i};
  end
  if model_external(i)
    report.external_data_ids{end+1, 1}  = data_ids{j};
    report.external_model_ids{end+1, 1} = model_ids{i};
  end
end

report.n_mapped     = length(report.mapped_data_ids);
report.has_problems = ~isempty(report.not_found) || ~isempty(report.ambiguous) || ...
                      ~isempty(report.compartment_ignored) || ~isempty(report.duplicate);
