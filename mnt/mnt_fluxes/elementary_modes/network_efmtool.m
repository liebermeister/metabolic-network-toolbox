function [C, mnet] = network_efmtool(network, method, ind_ignore, zv)

% [C,mnet] = network_efmtool(network, method, ind_ignore, zv)
%
% Compute the elementary flux modes by call the efmtool function "CalculateFluxModes"
%
% Function arguments:
%   method: 'internal': (default) consider only internal metabolites (-> all stationary modes)
%           'total': consider also external metabolites (-> all stationary, eba-unfeasible modes)
%   ind_ignore: indices of reactions (e.g. biomass production) to be neglected
%   zv (optional [used if non-empty]): vector of flux benefit values.
%      if zv is given, only modes without a flux benefit are given back

a = pwd;

if ~exist('CalculateFluxModes','file')
  error('Please install efmtool for computing elementary flux modes');
end
  
cd /home/wolfram/matlab_fixed/packages/efmtool/

eval(default('method','''internal''','ind_ignore','[]','zv','[]'));

switch method,
  case 'internal',
    ind_balances = find(network.external==0);
  case 'total',
    ind_balances = 1:length(network.metabolites);
  otherwise, 
    error('unknown method');
end

ind_reactions = setdiff(1:length(network.actions),ind_ignore);

stru.stoich          = full(network.N(ind_balances,ind_reactions));
stru.reversibilities = full(network.reversible(ind_reactions));

use_KEGG_id = 0;
if isfield( network,'metabolite_KEGGID'),
  if length(unique(network.metabolite_KEGGID)) == length(network.metabolite_KEGGID),
    use_KEGG_id = 1;
  end
end

if use_KEGG_id,
   stru.metaboliteNames = network.metabolite_KEGGID(ind_balances);
else,
   stru.metaboliteNames = network.metabolites(ind_balances);
end

%stru.reactionNames   = network.reaction_KEGGID(ind_reactions);
n_max = length(ind_reactions); list=cell(n_max,1); for it=1:n_max,  list{it} = ['R' num2str(it)]; end
stru.reactionNames = list;

ind_keep = find(sum(abs(stru.stoich),2)>0);
stru.stoich = stru.stoich(ind_keep,:);
stru.metaboliteNames = stru.metaboliteNames(ind_keep);

if ~isempty(zv),
  stru.stoich          = full([stru.stoich; column(zv)']);
  stru.metaboliteNames = [stru.metaboliteNames; {'Flux_benefit'}]; 
  if size(zv,2)>1, stru.metaboliteNames = [stru.metaboliteNames; repmat({'Constraint'},size(zv,2)-1,1)]; end
end

opts = CreateFluxModeOpts('arithmetic', 'fractional');
mnet = CalculateFluxModes(stru,opts);

cd(a);

if ~all(all(abs(mnet.stoich*mnet.efms) < 1e-10)),
  error('Computed elementary modes violate stationarity condition');
end

C = zeros(length(network.actions),size(mnet.efms,2));
C(ind_reactions,:) = mnet.efms;

switch method,
  case 'internal',
    if norm(network.N(find(network.external==0),:) * C)>10^-10,
      error('Computed cycles are not stationary');
    end
  case 'total',
    if norm(network.N * C)>10^-10,
      error('Computed cycles are not stationary');
    end
end
