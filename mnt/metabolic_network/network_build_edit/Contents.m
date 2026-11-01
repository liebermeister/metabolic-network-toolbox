% Build and edit metabolic networks
%
%   network_add_empty
%   network_add_metabolites              - Add metabolites to a network
%   network_add_reaction                 - Add a reaction to a network
%   network_add_transport
%   network_build_from_sum_formulae      - Build matlab network structure from reactions contained in SBtab file
%   network_build_from_sum_formulae_list - Build matlab network structure from reaction list
%   network_check_regulation_matrix
%   network_choose                       - Choose a subnetwork, including all metabolites of its reactions
%   network_choose_from_list             - Define a subnetwork by interactively choosing metabolites from a list
%   network_construct                    - Build a network structure from a stoichiometric matrix
%   network_listchoose_actions
%   network_listchoose_metabolites
%   network_merge_metabolites            - Merge each metabolite in m2_list into the corresponding metabolite in m1_list
%   network_print_formulae               - Convert a network structure into reaction formulae
%   network_reduce_to_active_subnetwork  - Reduce a network to the reactions that carry flux
%   network_remove                       - Remove metabolites from a network
%   network_remove_action                - Remove reactions from a network
%   network_set_all_external
%   network_set_external                 - Set metabolites to be external
%   network_subnetwork                   - Choose a subnetwork with given metabolites and reactions
%   network_subnetwork_balanced          - Reduce a network to reactions with considerable flux, keeping it balanced
