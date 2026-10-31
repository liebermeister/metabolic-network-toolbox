% METABOLIC_NETWORK Data structure 'network' for metabolic network models
%
% A network is a Matlab struct. To construct one, use network_construct, or load a model
% with sbtab_to_network, network_sbml_import, or cobra_to_network.
% With m metabolites and n reactions, the fields are:
%
% Mandatory fields:
%   N            (m x n matrix)     stoichiometric matrix
%   metabolites  (m x 1 cell)       metabolite names
%   actions      (n x 1 cell)       reaction names
%   reversible   (n x 1 vector)     1 for reversible reactions, 0 otherwise
%   external     (m x 1 vector)     1 for external (fixed) metabolites, 0 for internal ones
%
% Optional fields:
%   regulation_matrix  (n x m matrix)  allosteric regulation of reactions by metabolites
%                                      (positive: activation, negative: inhibition)
%   kinetics           (struct)        rate laws, see 'help kinetics'
%   graphics_par       (struct)        layout and display settings, see below
%   name               (string)        model name
%   metabolite_names   (m x 1 cell)    long names
%   action_names       (n x 1 cell)    long names
%   formulae           (n x 1 cell)    reaction formulae
%   EC                 (n x 1 cell)    EC numbers
%   metabolite_KEGGID  (m x 1 cell)    KEGG compound identifiers
%   reaction_KEGGID    (n x 1 cell)    KEGG reaction identifiers
%   metabolite_id      (m x 1 cell)    unique identifiers of metabolites
%   action_id          (n x 1 cell)    unique identifiers of reactions
%   parameter_id       (cell)          unique identifiers of parameters
%   K                  (matrix)        kernel matrix
%   L                  (matrix)        link matrix
%   NR                 (matrix)        stoichiometric matrix of independent metabolites
%   s_init             (m x 1 vector)  initial concentrations
%
% Field graphics_par (created by netgraph_make_graph). Main fields:
%   x                  (2 x (m+n) matrix)  positions of metabolites, followed by reactions
%   fixed              (1 x (m+n) vector)  1 for positions that are not moved by relaxation
%   metnames, actnames (cell)              names shown in graphics
%   For the other display settings, see 'help netgraph_draw'.
%
% Other data structures:
%   kinetics      rate laws and their parameters ('help kinetics')
%   kinetic_data  measured or estimated parameter values ('help kinetic_data')
%   state_data    metabolite levels, fluxes, and enzyme levels ('help state_data')
%
% Functions: see the folders network_build_edit, network_inspect, and network_format_conversion
% To check a network structure, use network_check and network_check_consistency.
