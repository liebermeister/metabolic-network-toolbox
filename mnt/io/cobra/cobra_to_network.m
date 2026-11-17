function [network, ind_keepm, ind_keepr] = cobra_to_network(model,options)

% COBRA_TO_NETWORK Convert a COBRA model to a network structure
%
% [network, ind_keepm, ind_keepr] = cobra_to_network(model,options)
%
% model    COBRA model struct (fields mets, rxns, S, metNames, rxnNames; optional rev, lb, c, subSystems,
%          metKeggId), or the name of a .mat file with the variable 'model'
% options  (optional) struct with fields (defaults in brackets):
%
%   Exchange reactions and external metabolites
%   .flag_restrict            [0]  remove the metabolites in .rem and the reactions of subsystem 'Exchange'
%   .rem                      [{}] metabolites removed if flag_restrict = 1
%   .remove_reactions_prefix  ['']  remove reactions whose IDs start with this string (e.g. 'EX_' for
%                                   exchange reactions)
%   .external_method          ['exchange'] how external metabolites are determined:
%                                   'exchange'     metabolites in reactions of subsystem 'Exchange'
%                                   'compartment'  metabolites in the compartments .external_compartments
%                                                  (ID suffix '[e]' etc.)
%   .external_compartments    [{'e'}] compartment labels for external_method = 'compartment'
%   .external_single_reaction [0]  also set metabolites external that occur in only one (remaining) reaction
%   .flag_set_external        [0]  network_set_external(network,1) at the end
%
%   Identifiers and names
%   .metabolite_id_replacements [{}] cell array {old1, new1; old2, new2; ...}: substrings replaced in all
%                                   metabolite IDs (strrep, in this order), e.g. {'-L',''; '-D',''} to remove
%                                   hyphens (which cause problems in SBtab reaction formulae)
%   .replace_metabolites      [{}] metabolite IDs (after the substring replacements) replaced by ...
%   .replace_metabolites_by   [{}] ... these IDs
%   .add_kegg_ids             [0]  add KEGG IDs (fields metabolite_KEGGID, reaction_KEGGID) from the BiGG IDs
%                                   (convert_bigg_to_kegg)
%   .id_mapping_file          [[]] mapping file for convert_bigg_to_kegg (default: the MNT mapping table)
%   .bigg_ids                 [0]  BiGG IDs = model IDs (fields metabolite_BIGGID, reaction_BIGGID)
%   .names_for_plots          [0]  fields metabolite_NameForPlots, reaction_NameForPlots = IDs
%   .reversible_and_formulae  [1]  fields reversible (from model.rev or model.lb) and formulae
%   .full_matrix              [0]  stoichiometric and regulation matrix as full (not sparse) matrices
%
% Biomass reactions: the reactions of the COBRA objective (model.c ~= 0) are stored in
% network.biomass_reactions (SBtab: Reaction table attribute BiomassReactionID), if model.c exists.
%
% Since 2026-11-17, this function also covers the conversion previously done by the IMF function
% imf_convert_cobra_model (E. coli core model in IMF):
%   cobra_to_network(model, struct('remove_reactions_prefix','EX_', 'external_method','compartment', ...
%     'external_single_reaction',1, 'metabolite_id_replacements',{{'-C','_C'; '-L',''; '-D',''}}, ...
%     'add_kegg_ids',1, 'bigg_ids',1, 'names_for_plots',1, 'reversible_and_formulae',0, 'full_matrix',1))
%
% See also: network_to_cobra, cobra_to_sbtab

options_default = struct('flag_restrict',0,'flag_set_external',0,'remove_reactions_prefix','', ...
                         'external_method','exchange','external_single_reaction',0, ...
                         'add_kegg_ids',0,'id_mapping_file',[],'bigg_ids',0,'names_for_plots',0, ...
                         'reversible_and_formulae',1,'full_matrix',0);
options_default.rem                        = {};
options_default.external_compartments      = {'e'};
options_default.metabolite_id_replacements = {};
options_default.replace_metabolites        = {};
options_default.replace_metabolites_by     = {};

eval(default('options','struct'));
options = join_struct(options_default,options);

if ischar(model)
  a = load(model);
  model = a.model;
end

% ------------------------------------------------------------
% metabolites and reactions to keep

if isfield(model,'subSystems')
  ind_exchange = find(strcmp(model.subSystems,'Exchange'));
else
  ind_exchange = [];
end

if options.flag_restrict
  display('Removing some metabolites and reactions');
  ind_rem   = label_names(options.rem,model.mets);
  ind_keepm = setdiff(1:length(model.mets),ind_rem);
  ind_keepr = find(~strcmp(model.subSystems,'Exchange'));
else
  ind_keepm = 1:length(model.mets);
  ind_keepr = 1:length(model.rxns);
end
ind_keepm = column(ind_keepm);
ind_keepr = column(ind_keepr);

if length(options.remove_reactions_prefix)
  ind_keepr = ind_keepr(~startsWith(model.rxns(ind_keepr), options.remove_reactions_prefix));
end

% ------------------------------------------------------------
% network

network.metabolites = column(model.mets(ind_keepm));
network.actions     = column(model.rxns(ind_keepr));
network.N           = model.S(ind_keepm,ind_keepr);
if options.full_matrix
  network.N = full(network.N);
end

switch options.external_method
  case 'exchange'
    external = zeros(length(model.mets),1);
    external(find(sum(abs(model.S(:,ind_exchange)),2))) = 1;
    network.external = external(ind_keepm);
  case 'compartment'
    network.external = zeros(length(network.metabolites),1);
    for it = 1:length(options.external_compartments)
      network.external(endsWith(network.metabolites, ['[' options.external_compartments{it} ']'])) = 1;
    end
  otherwise
    error('Unknown external_method ''%s''', options.external_method);
end

if options.external_single_reaction
  network.external(sum(network.N~=0,2) == 1) = 1;
end

if options.full_matrix
  network.regulation_matrix = zeros(size(network.N'));
else
  network.regulation_matrix = sparse(zeros(size(network.N')));
end
if options.reversible_and_formulae
  if isfield(model,'rev')
    network.reversible = model.rev(ind_keepr);
  elseif isfield(model,'lb')
    network.reversible = double(model.lb(ind_keepr) < 0);
  end
  network.formulae = network_print_formulae(network);
end
network.metabolite_names = column(model.metNames(ind_keepm));
network.reaction_names   = column(model.rxnNames(ind_keepr));
if options.names_for_plots
  network.metabolite_NameForPlots = network.metabolites;
  network.reaction_NameForPlots   = network.actions;
end
if isfield(model,'metKeggId')
  network.metabolite_KEGGID = column(model.metKeggId(ind_keepm));
end

% biomass reactions: reactions in the objective function of the COBRA model
if isfield(model, 'c')
  network.biomass_reactions = intersect(model.rxns(model.c(:) ~= 0), network.actions, 'stable');
else
  network.biomass_reactions = {};
end

% ------------------------------------------------------------
% metabolite IDs

for it = 1:size(options.metabolite_id_replacements,1)
  network.metabolites = strrep(network.metabolites, options.metabolite_id_replacements{it,1}, options.metabolite_id_replacements{it,2});
end

if options.add_kegg_ids
  [network.metabolite_KEGGID, network.reaction_KEGGID] = convert_bigg_to_kegg(network.metabolites, network.actions, options.id_mapping_file);
end

for it = 1:length(options.replace_metabolites)
  network.metabolites{network_find_metabolites(network,options.replace_metabolites{it})} = options.replace_metabolites_by{it};
end

if options.bigg_ids
  network.metabolite_BIGGID = network.metabolites;
  network.reaction_BIGGID   = network.actions;
end

if options.flag_set_external
  network = network_set_external(network,1);
end
