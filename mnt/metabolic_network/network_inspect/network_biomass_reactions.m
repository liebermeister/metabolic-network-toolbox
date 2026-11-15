function biomass_reactions = network_biomass_reactions(network, default_reactions)

% NETWORK_BIOMASS_REACTIONS Biomass reactions of a network (reaction IDs)
%
% biomass_reactions = network_biomass_reactions(network, default_reactions)
%
% network            MNT network struct
% default_reactions  (optional) cell array of reaction IDs, used if network.biomass_reactions
%                    does not exist (default {})
%
% biomass_reactions  column cell array of reaction IDs (network.actions): network.biomass_reactions
%                    if this field exists (even if empty), otherwise those of default_reactions
%                    that occur in the network
%
% The field network.biomass_reactions (cell array of reaction IDs) is read from and written to the
% attribute BiomassReactionID of the SBtab Reaction table (comma-separated list), see sbtab_to_network
% and network_to_sbtab.

if ~exist('default_reactions', 'var'), default_reactions = {}; end
if ischar(default_reactions), default_reactions = {default_reactions}; end

if isfield(network, 'biomass_reactions')
  biomass_reactions = cellstr(network.biomass_reactions);
else
  biomass_reactions = default_reactions(ismember(default_reactions, network.actions));
end
biomass_reactions = biomass_reactions(:);
