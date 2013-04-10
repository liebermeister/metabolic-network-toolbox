function [subnetwork,indm,indr] = network_choose(network,indm,indr)

% [subnetwork,indm,indr] = network_choose(network,indm,indr)
% 
% Choose subnetwork containing metabolites/reactions with
% the (old) indices indm (for metabolites) and indr (for reactions) 
% and all other metabolites taking part in these reactions

if ~exist('indr','var'), indr = 'add_double'; end

if isstr(indr),
  switch indr,
    case 'add_all',
      disp('Adding all reactions connected to any of the metabolites');
      indr = find(sum(abs(network.N(indm,:)),1)); 
    case 'add_double',
      disp('Adding all reactions connected to at least two metabolites');
      indr = find(sum((network.N(indm,:)~=0),1)>1);
  end
end

indm = union(indm,find(sum(abs(network.N(:,indr)),2)));   

dummi                  = zeros(1,length(network.metabolites));
subnetwork.metabolites = network.metabolites(indm);
subnetwork.actions     = network.actions(indr);
subnetwork.N           = network.N(indm,indr);
subnetwork.reversible  = network.reversible(indr);
dummi(find(network.external))=1;
subnetwork.external   = bit_vector(find(dummi(indm)),length(subnetwork.metabolites));

if isfield(network,'metabolite_KEGGID'), subnetwork.metabolite_KEGGID = network.metabolite_KEGGID(indm); end
if isfield(network,'reaction_KEGGID'),   subnetwork.reaction_KEGGID   = network.reaction_KEGGID(indr); end
if isfield(network,'EC'),                subnetwork.EC                = network.EC(indr); end
if isfield(network,'formulae'),          subnetwork.formulae          = network.formulae(indr); end
if isfield(network,'index'),             subnetwork.index             = network.index(indr); end

if isfield(network,'regulation_matrix'),
  subnetwork.regulation_matrix = network.regulation_matrix(indr,indm);
end

if isfield(network,'graphics_par'), 
  x = network.graphics_par.x;
  x_new = [x(:,indm), x(:,length(network.metabolites)+indr)];
  subnetwork=netgraph_make_graph(subnetwork); 
  subnetwork.graphics_par.x = x_new;
end

if isfield(network,'kinetics'),

  subnetwork.kinetics    = network.kinetics;
 
  switch subnetwork.kinetics.type,
    case 'mass-action',
      subnetwork.kinetics.k_fwd = subnetwork.kinetics.k_fwd(indr);
      subnetwork.kinetics.k_bwd = subnetwork.kinetics.k_bwd(indr);
    case 'numeric',
      subnetwork.kinetics.use_only_act = indr; 
      subnetwork.kinetics.use_only_met = indm; 
      subnetwork.kinetics.n_met_tot    = length(indm);
    case 'saturated',
      subnetwork.kinetics.reactions = subnetwork.kinetics.reactions(indr); 
    case 'numeric',
      subnetwork.kinetics = submodel_numeric(subnetwork.kinetics,length(subnetwork.metabolites),indices_met_sub,indr);
    case {'cs','ms','rp','ma','fm'},
      subnetwork.kinetics = modular_reduce_to_subnetwork(network.kinetics,indm,indr);
    otherwise,
      error(sprintf('Kinetics type %s not supported in subnetwork selection',subnetwork.kinetics.type));
  end
end
