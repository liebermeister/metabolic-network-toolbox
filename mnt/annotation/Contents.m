% Annotation: mapping between identifiers of metabolites and reactions
%
%   convert_bigg_to_kegg            - Convert BiGG ids of metabolites and reactions to KEGG ids
%   convert_kegg_to_bigg            - Convert KEGG ids of metabolites and reactions to BiGG ids (with optional extra mappings)
%   bigg_kegg_mapping_load          - Load a BiGG-KEGG mapping table
%   bigg_kegg_strip_metabolite_ids  - Remove compartment suffixes from BiGG metabolite ids
%   bigg_metabolite_id_base         - Split BiGG-style metabolite ids into base id and compartment
%   bigg_match_metabolite_ids       - Match metabolite ids of a data set to those of a model (with report)
%   bigg_match_metabolite_ids_display - Print the report of bigg_match_metabolite_ids
%
% Data
%   bigg-kegg-mapping.tsv           - Mapping between BiGG and KEGG identifiers (BiGG reactions with their
%                                     metabolites, stoichiometries, KEGG reaction and compound ids, EC numbers);
%                                     from the IMF toolbox (resources/data/annotation-data/BiggMapping.tsv, 2023)
