function mnt_save_model_and_state(model_outfile, network, v, c, keq, networkG, options)

% MNT_SAVE_MODEL_AND_STATE Save a model with its reference state (or a network without state) in .mat, SBtab and SBML files
%
% mnt_save_model_and_state(model_outfile, network, v, c, keq, networkG, options)
%
% Counterpart of mnt_load_model_and_state. Writes the files
%   [model_outfile '.mat']  Matlab file with the variables network, networkG, v, c, keq
%                           (and further variables given in options.mat_variables, e.g. A, dG0)
%   [model_outfile '.tsv']  SBtab file (network_to_sbtab: tables Reaction, Compound, Compartment,
%                           concentrations, fluxes, reaction Gibbs free energies, equilibrium constants),
%                           with the table 'NetworkLayout' for network graphics (netgraph_layout_table)
%                           if networkG is given
%   [model_outfile '.xml']  SBML file (sbml_network_state_save: network, concentrations, fluxes,
%                           equilibrium constants, reaction Gibbs free energies; no layout)
%
% model_outfile  file name without extension (e.g. '.../glucose-aerobic/Model'); the folder must exist
% network        network struct; kinetics (network.kinetics) are not written; for a network without state
%                and with networkG, the layout of the full network (network.graphics_par) is not written either
% v, c, keq      fluxes, concentrations, equilibrium constants; all three empty ([]): network without
%                state (e.g. a basic network model; then the .mat file contains empty vectors)
% networkG       (optional) network with layout for graphics (e.g. from netgraph_make_layout_network); [] for none
% options        (optional) struct with fields
%   .mat_variables  struct with further variables for the .mat file          (default: none)
%   .formats        files to write, subset of {'mat','tsv','xml'}             (default: all three)
%   .model_name     model name in the SBML file (default: name of the folder that contains the file;
%                   hyphens replaced by underscores)
%
% mnt_load_model_and_state reads the files; for SBtab files, networkG is rebuilt from the table
% NetworkLayout (netgraph_layout_from_table).
%
% (Until 2026-11-17 the IMF functions imf_save_model_and_state and imf_save_network_model.)
%
% See also: mnt_load_model_and_state, network_to_sbtab, sbml_network_state_save, netgraph_layout_table

eval(default('networkG','[]','options','struct'));
[~, folder_name] = fileparts(fileparts(model_outfile));
options = join_struct(struct('mat_variables', struct, 'formats', {{'mat','tsv','xml'}}, 'model_name', folder_name), options);

if isfield(network, 'kinetics'), network = rmfield(network, 'kinetics'); end

has_state = ~isempty(v) || ~isempty(c) || ~isempty(keq);
% network without state: layout of the full network (graphics_par) not saved if networkG is given
% (the layout is in the table NetworkLayout)
if ~has_state && ~isempty(networkG) && isfield(network, 'graphics_par'), network = rmfield(network, 'graphics_par'); end
if has_state && ~isempty(keq) && ~isempty(c)
  dmu = -RT * (log(keq) - network.N' * log(c));
else
  dmu = [];
end

% .mat file
if ismember('mat', options.formats)
  S = options.mat_variables;
  S.network  = network;
  S.networkG = networkG;
  S.v        = v;
  S.c        = c;
  S.keq      = keq;
  save([model_outfile '.mat'], '-struct', 'S');
end

% SBtab file
if ismember('tsv', options.formats)
  sbtab = network_to_sbtab(network, struct('c', c, 'v', v, 'dmu', dmu, 'keq', keq, 'verbose', 0));
  if ~isempty(networkG)
    sbtab = sbtab_document_add_table(sbtab, 'NetworkLayout', netgraph_layout_table(networkG));
  end
  sbtab_document_save(sbtab, [model_outfile '.tsv'], 0, 0);
end

% SBML file
if ismember('xml', options.formats)
  sbml_network_state_save(network, c, v, [model_outfile '.xml'], dmu, struct('name', strrep(options.model_name,'-','_'), 'verbose', 0), keq);
end
