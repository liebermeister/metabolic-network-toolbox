%[network,s] = network_sbml_import(s,dirname,verbose)
%
%Import a network from SBML format.
%
%s can be:
%- a SBML structure (as defined by the SBMLToolbox)
%- a filename. In this case, the firectory name has to be given 
%  in the argument 'dirname'
%
%Kinetic laws are represented in kinetics of the 'numeric'
%type (see 'network_velocities')
%
%This function requiers the SBMLToolbox to be installed
%(see http://sbml.org/software/sbmltoolbox/)

function [network,s] = network_sbml_import(s,dirname,verbose)

if ~exist('TranslateSBML','file'),
  error('Please install the SBML Toolbox (http://sbml.org/Software/SBMLToolbox) - Otherwise the  SBML import/export functions do not work.');
end

eval(default('verbose','0'));

create_annotations = 0;   % CONSTRUCT UNIQUE ANNOTATIONS ?
graphics_flag      = 0;   % SHOW GRAPHICS?
import_names       = 1;   % import_names from model
import_compartments= 0;   % import_names from model

if isstr(s),
  if exist('dirname','var'); cd(dirname); end 
  s = TranslateSBML(s);
end

% -----------------------------------------------------------------------------
% --- analyse list of species -> 

%  lists 'metabolites', 'metabolite_names', 
% 'metabolite_compartments', 's_init', 'a_init', 'external', 'metabolite_units'

clear metabolites 
for it=1:length(s.species),
  metabolites{it,:} = s.species(it).id;
  metabolite_names{it,:} = s.species(it).name;
  metabolite_compartments{it,:} = s.species(it).compartment;
  external(it,:) = s.species(it).boundaryCondition;
  s_init(it,:)   = s.species(it).initialConcentration;
  a_init(it,:)   = s.species(it).initialAmount;
  ind = findstr(s.species(it).annotation,'kegg.compound'); 
  if ind, metabolite_kegg_id{it,:} = s.species(it).annotation(ind+14:ind+19); 
  else,
    ind = findstr(s.species(it).annotation,'kegg'); 
    if ind, metabolite_kegg_id{it,:} = s.species(it).annotation(ind+5:ind+10); 
    else,
      metabolite_kegg_id{it,:} = '';
    end
  end
end

% -----------------------------------------------------------------------------
% --- analyse list of reactions

actions    = [];
reversible = [];

for it=1:length(s.reaction),
  reaction_id{it,:} = s.reaction(it).id; 
  reaction_name{it,:} = s.reaction(it).name; 
  reaction_fast{it,:} = s.reaction(it).fast;
  reversible(it,:) = s.reaction(it).reversible;
  reaction_annotation{it,:} = s.reaction(it).annotation;
  if isempty(s.reaction(it).kineticLaw),
    reaction_formula{it,:}    = [];
    reaction_math{it,:}       = [];
    reaction_parameter{it,:}  = [];
  else,
    reaction_formula{it,:}    = s.reaction(it).kineticLaw.formula;
    reaction_math{it,:}       = s.reaction(it).kineticLaw.math;
    reaction_parameter{it,:}  = s.reaction(it).kineticLaw.parameter;
  end
  listOfReactants{it,:} = s.reaction(it).reactant;
  listOfProducts{it,:}  = s.reaction(it).product;
  listOfModifiers{it,:} = s.reaction(it).modifier;
  ind = findstr(s.reaction(it).annotation,'kegg.reaction'); 
  if ind, 
    reaction_kegg_id{it,:} = s.reaction(it).annotation(ind+14:ind+19);
  else,
    reaction_kegg_id{it,:} = '';
  end

end

actions = reaction_id;
if isempty(reversible), reversible = ones(size(actions)); else, reversible=reversible; end

% -----------------------------------------------------------------------------
% analyse kinetics -> lists stoichiometries, pstoichiometries


for it=1:length(listOfReactants),
  substrate{it}={};
  sstoichiometries{it}=[];
  for itt = 1:length(listOfReactants{it});
    substrate{it}= [substrate{it} {listOfReactants{it}(itt).species}];
    sstoichiometries{it}= [sstoichiometries{it} listOfReactants{it}(itt).stoichiometry];
  end
end
  
for it=1:length(listOfProducts),
  product{it}={};
  pstoichiometries{it}=[];
  for itt = 1:length(listOfProducts{it});
    product{it}= [product{it} {listOfProducts{it}(itt).species}];
    pstoichiometries{it}= [pstoichiometries{it} listOfProducts{it}(itt).stoichiometry];
  end
end

for it=1:length(listOfModifiers),
  modifier{it}={};
  for itt = 1:length(listOfModifiers{it});
    modifier{it}= [modifier{it} {listOfModifiers{it}(itt).species}];
  end
end

% ----------------------------------------------------------------

N=zeros(length(metabolites),length(reaction_id));

for it=1:length(s.reaction),
  substrate_indices =  label_names(substrate{it},metabolites,'single');
  if length(substrate_indices), N(substrate_indices,it) = N(substrate_indices,it)-sstoichiometries{it}'; end 
  product_indices =  label_names(product{it},metabolites,'single');
  if length(product_indices),
    N(product_indices,it) = N(product_indices,it) + pstoichiometries{it}' ; end
end

% --------------------------------------------------------------------------------------------

network.metabolites  = metabolites;
network.actions      = reaction_id;
network.external     = double(external);
network.reversible   = reversible;
network.N            = N;

network.metabolite_KEGGID = metabolite_kegg_id;
network.reaction_KEGGID   = reaction_kegg_id;

if import_names,

if exist('metabolite_names','var'), network.metabolite_names  = metabolite_names; 
else, network.metabolite_names  = network.metabolites; end

if exist('reaction_name','var'), network.reaction_names  = reaction_name; 
else, network.reaction_names  = network.actions; end

end

if import_names,
  if exist('metabolite_names','var'), network.metabolite_compartments  = metabolite_compartments; 
  else, network.metabolite_compartments  = repmat({''},length(network.metabolites),1); end
end

% ------------------------------------------------------------
% set reaction kinetics

network.regulation_matrix = zeros(size(network.N))';
  
for it1 = 1:length(modifier),
  for it2=1:length(modifier{it1}),
    mm = label_names(modifier{it1}(it2),metabolites);
    if mm,
      network.regulation_matrix(it1,mm)=1;
    end
  end
end

if ~isempty(reaction_formula{1}),
  network.kinetics.type = 'numeric';
  network.kinetics.reactions = {};
  all_parameters = {};
  for it = 1:length(actions),
    network.kinetics.reactions{it}.string = strrep(reaction_math{it},'pow','power');
    for itt = 1:length(reaction_parameter{it}),
      this_parameter = reaction_parameter{it}(itt);
      network.kinetics.reactions{it}.parameters{itt}.name = this_parameter.id;
      network.kinetics.reactions{it}.parameters{itt}.value = this_parameter.value;
    end
    if length(reaction_parameter{it}),
      network.kinetics.reactions{it}.parameters = network.kinetics.reactions{it}.parameters';
      else,
      network.kinetics.reactions{it}.parameters = {};
    end
    end
  network.kinetics.reactions = network.kinetics.reactions';
end


if ~isempty(reaction_formula{1}),
  z=1;
  network.kinetics.parameter_values = [];
  network.kinetics.parameters = [];

  for it=1:length(network.actions),
    par = network.kinetics.reactions{it}.parameters;
    for it2 = 1: length(par),
      network.kinetics.parameters{z}   = ['v' num2str(it) '_' par{it2}.name];
      network.kinetics.parameter_values(z)   = par{it2}.value;
      network.kinetics.reactions{it}.parameters{it2}.index = z;
      network.kinetics.reactions{it}.parameters{it2} = rmfield(network.kinetics.reactions{it}.parameters{it2},'value');
      z=z+1;
    end
  end
  network.kinetics.parameter_values = network.kinetics.parameter_values';
  network.kinetics.parameters = network.kinetics.parameters';
elseif verbose, 
  fprintf('Warning: no kinetics found in SBML file.\n');
end 

if exist('s_init','var'), network.s_init       = s_init; 
else,
if exist('a_init','var'), network.s_init       = a_init; end
end

if graphics_flag,
  if length(network.actions),
    netgraph_draw(network);
  end
end


% --------------------------------------------------------------

%if create_annotations,

%for it=1:length(network.metabolites),
%  network.metabolite_annotation{it}.object      = network.metabolites{it};  
%  network.metabolite_annotation{it}.is_reaction = 0;
%  network.metabolite_annotation{it}.quantity    = 'concentration';
%  network.metabolite_annotation{it}.unit        = 'mM';
%end

%for it=1:length(network.actions),
%  network.action_annotation{it}.object        = network.actions{it};
%  network.action_annotation{it}.is_reaction   =  1;
%  network.action_annotation{it}.quantity      = 'flux';
%  network.action_annotation{it}.unit          = 'mM/min';
%end

%if ~isempty(reaction_formula),
%  z=1;
%  for it=1:length(network.actions),
%    par = network.kinetics.reactions{it}.parameters;
%    for it2 = 1: length(par),
%    network.parameter_annotation{z}.quantity      = 'reaction parameter';
%    network.parameter_annotation{it}.is_reaction  = 0;
%    network.parameter_annotation{it}.unit         = 'unknown';
%    network.parameter_annotation{z}.object        = par{it2}.name;
%    network.parameter_annotation{z}.value         = par{it2}.value;
%    network.parameter_annotation{z}.appears_in_reaction = it;
%    z=z+1;
%    end
%  end
%end
