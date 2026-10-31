% NETWORK_SBML_EXPORT Export a network to SBML
%
% SBMLmodel = network_sbml_export(network, verbose, name, filename, notes, sbml_level, sbml_version)
%
% Export network (data type from mnt) to libSBML's SBMLModel data type
%
% ARGUMENTS
%  network   (network structure) network
%  name      (string, optional)  model name
%  notes     (string, optional)
%
% If the argument filename (string) is given, the model is directly written
% to an SBML (.xml) file (in the present directory; no path can be given)
%
% You can also use OutputSBML(SBMLmodel,filename) to export the model
% to an SBML (.xml) file
%
% This function requires the Matlab interface of libSBML
% (TranslateSBML and OutputSBML; see https://sbml.org/software/libsbml/)
%
% See also: network_sbml_import, network_to_sbtab

function SBMLmodel = network_sbml_export(network,verbose,name,filename,notes,sbml_level,sbml_version)

if ~exist('TranslateSBML','file')
  error('Please install the Matlab interface of libSBML (https://sbml.org/software/libsbml/) - Otherwise the SBML import/export functions do not work.');
end

eval(default('verbose','0','name','''Model''','notes','''''','sbml_level','2','sbml_version','2','filename','[]'));

if length(filename)==0, filename = []; end

name = strrep(name,' ','_');
name = strrep(name,'.','_');


% -------------------------------------------------------------
% Adjust names to SBML conventions

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

if length(network.metabolites) ~= length(unique(network.metabolites))
  warning('Error in conversion of metabolite names');
  for it=1:length(network.metabolites)
    network.metabolites{it} = ['M' num2str(it) '_' network.metabolites{it}];
  end
end

if length(unique(network.actions))<length(network.actions)
  for it=1:length(network.actions)
    network.actions{it} = ['A' num2str(it) '_' network.actions{it}];
  end
end


% -------------------------------------------------------------
% SBML object
%
% The model is built as a Matlab struct in the format used by libSBML's
% TranslateSBML and OutputSBML. Empty structs of the right format are
% obtained from a minimal SBML file (see sbml_templates below).

t = sbml_templates(sbml_level,sbml_version);

SBMLmodel      = t.model;
SBMLmodel.name = name;
SBMLmodel.id   = name;

compartment      = t.compartment;
compartment.id   = 'compartment';
compartment.name      = 'compartment';
compartment.size      = 1;
compartment.isSetSize = 1;
SBMLmodel.compartment = compartment;

% notes
% Model %s written by the Metabolic Network Toolbox</p>\n',name)];


% ----------------------------------------------------------------
% metabolites

if verbose, fprintf('  Converting the metabolites: '); end

species = t.species([]);

for it = 1:length(network.metabolites)
  if verbose,  fprintf('%d ',it); end
  species(it)                   = t.species;
  species(it).id                = network.metabolites{it};
  species(it).name              = network.metabolites{it};
  species(it).compartment       = 'compartment';
  species(it).boundaryCondition = double(network.external(it));
  if isfield(network,'s_init')
    species(it) = set_initial_concentration(species(it),network.s_init(it));
  else
    species(it) = set_initial_concentration(species(it),0);
  end
end


% ----------------------------------------------------------------
% kinetic laws

if isfield(network,'kinetics')
  switch network.kinetics.type

    case {'ms','cs','ds','rp'}

      if ~isfield(network.kinetics,'c')
        network.kinetics.c = nan * ones(size(network.metabolites));
      end
      for it = 1:length(network.metabolites)
        species(it) = set_initial_concentration(species(it),network.kinetics.c(it));
      end

      metnames = network.metabolites;
      kk       = network.kinetics;

      if verbose, fprintf('\n  Converting the rate laws: '); end

      for it = 1:length(network.actions)

        if verbose, fprintf('%d ',it); end

        % THIS IS A COPY FROM network_get_kinetics_strings
        r_name = ['R' num2str(it) ];

        sub = find(network.N(:,it) < 0);
        pro = find(network.N(:,it) > 0);
        rea = find(kk.KM(it,:)~=0);
        act = find(kk.KA(it,:)~=0);
        inh = find(kk.KI(it,:)~=0);

        m_sub           = abs(network.N(sub,it));
        m_pro           = abs(network.N(pro,it));

        switch network.kinetics.type
          case 'ms'
            formula    = ms_get_formula(r_name,metnames,sub,pro,act,inh,m_sub,m_pro);
            kk.KVratio = ms_KM_KVratio_to_Keq(network.N,kk.KM,kk.Keq);
          case 'cs'
            formula    = cs_get_formula(r_name,metnames,sub,pro,act,inh,m_sub,m_pro);
            kk.KVratio = ms_KM_KVratio_to_Keq(network.N,kk.KM,kk.Keq);
          case 'rp'
            formula    = ms_get_formula(r_name,metnames,sub,pro,act,inh,m_sub,m_pro);
            kk.KVratio = ms_KM_KVratio_to_Keq(network.N,kk.KM,kk.Keq);
          case 'ds'
            formula    = cs_get_formula(r_name,metnames,sub,pro,act,inh,m_sub,m_pro);
            kk.KVratio = ms_KM_KVratio_to_Keq(network.N,kk.KM,kk.Keq);
        end

        if ~isfield(kk,'u')
          kk.u = nan * ones(size(network.actions));
        end

        par_ids    = {['u_' r_name ], ['kC_' r_name ], ['kEQ_' r_name ]};
        par_values = [kk.u(it), kk.KV(it), kk.Keq(it)];
        for itt = 1:length(rea)
          par_ids{end+1}    = ['kM_' r_name '_' metnames{rea(itt)} ];
          par_values(end+1) = kk.KM(it,rea(itt));
        end
        for itt = 1:length(act)
          par_ids{end+1}    = ['kA_' r_name '_' metnames{act(itt)} ];
          par_values(end+1) = kk.KA(it,act(itt));
        end
        for itt = 1:length(inh)
          par_ids{end+1}    = ['kI_' r_name '_' metnames{inh(itt)} ];
          par_values(end+1) = kk.KI(it,inh(itt));
        end

        kinetic_law{it} = make_kinetic_law(t, formula, par_ids, par_values);
      end

    case 'kinetic_strings'

      for it = 1:length(network.actions)
        r = network.kinetics.reactions{it};
        kinetic_law{it} = make_kinetic_law(t, strrep(r.string,'power','pow'), r.parameters, r.parameter_values);
      end

    case 'numeric'

      if isfield(network.kinetics,'velocity_strings')
        my_strings = feval(network.kinetics.velocity_strings);
        for itt = 1:length(my_strings)
          kinetic_law{itt} = make_kinetic_law(t, my_strings{itt}, {}, []);
        end

        parnames = fieldnames(network.kinetics.parameters);
        for itt = 1:length(parnames)
          SBMLmodel.parameter(itt) = make_parameter(t.parameter, parnames{itt}, network.kinetics.parameters.(parnames{itt}));
        end
      end

    otherwise, disp('Kinetics export has not been implemented for this kinetics so far.');

  end
end

SBMLmodel.species = species;


% -----------------------------------------------------------
% reactions

if verbose, fprintf('\n  Converting the reactions: '); end

reactions = t.reaction([]);

for it = 1:length(network.actions)

  if verbose,  fprintf('%d ',it); end
  reaction            = t.reaction;
  reaction.id         = network.actions{it};
  reaction.name       = network.actions{it};
  reaction.reversible = double(network.reversible(it));
  if exist('kinetic_law','var')
    reaction.kineticLaw = kinetic_law{it};
  end

  sub_ind = find(network.N(:,it)<0);
  for itt = 1:length(sub_ind)
    reaction.reactant(itt)               = t.speciesReference;
    reaction.reactant(itt).species       = network.metabolites{sub_ind(itt)};
    reaction.reactant(itt).stoichiometry = abs(network.N(sub_ind(itt),it));
  end

  prod_ind = find(network.N(:,it)>0);
  for itt = 1:length(prod_ind)
    reaction.product(itt)               = t.speciesReference;
    reaction.product(itt).species       = network.metabolites{prod_ind(itt)};
    reaction.product(itt).stoichiometry = abs(network.N(prod_ind(itt),it));
  end

  if isfield(network,'regulation_matrix')
    mod_ind = find(network.regulation_matrix(it,:));
    for itt = 1:length(mod_ind)
      reaction.modifier(itt)         = t.modifier;
      reaction.modifier(itt).species = network.metabolites{mod_ind(itt)};
    end
  end

  reactions(it) = reaction;

end

SBMLmodel.reaction = reactions;

if verbose, fprintf('\n',it); end

if length(filename)
  OutputSBML(SBMLmodel,filename);
  fprintf('Wrote SBML file "%s"\n',filename);
end


% -------------------------------------------------------------------

function t = sbml_templates(sbml_level,sbml_version)

  % Empty libSBML structs (model, compartment, species, parameter, reaction,
  % species reference, modifier, kinetic law, and kinetic law parameter)
  % for a given SBML level and version, obtained by reading a minimal SBML
  % file with TranslateSBML

  if sbml_level < 3
    xml = sprintf(['<?xml version="1.0" encoding="UTF-8"?>\n' ...
                   '<sbml xmlns="http://www.sbml.org/sbml/level%d/version%d" level="%d" version="%d">\n' ...
                   '<model id="m">\n' ...
                   '<listOfCompartments><compartment id="c"/></listOfCompartments>\n' ...
                   '<listOfSpecies><species id="A" compartment="c"/></listOfSpecies>\n' ...
                   '<listOfParameters><parameter id="p"/></listOfParameters>\n' ...
                   '<listOfReactions><reaction id="r">\n' ...
                   '<listOfReactants><speciesReference species="A"/></listOfReactants>\n' ...
                   '<listOfModifiers><modifierSpeciesReference species="A"/></listOfModifiers>\n' ...
                   '<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>A</ci></math>\n' ...
                   '<listOfParameters><parameter id="k"/></listOfParameters></kineticLaw>\n' ...
                   '</reaction></listOfReactions>\n' ...
                   '</model>\n</sbml>\n'], sbml_level, sbml_version, sbml_level, sbml_version);
  else
    xml = sprintf(['<?xml version="1.0" encoding="UTF-8"?>\n' ...
                   '<sbml xmlns="http://www.sbml.org/sbml/level%d/version%d/core" level="%d" version="%d">\n' ...
                   '<model id="m">\n' ...
                   '<listOfCompartments><compartment id="c" constant="true"/></listOfCompartments>\n' ...
                   '<listOfSpecies><species id="A" compartment="c" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>\n' ...
                   '<listOfParameters><parameter id="p" constant="true"/></listOfParameters>\n' ...
                   '<listOfReactions><reaction id="r" reversible="true" fast="false">\n' ...
                   '<listOfReactants><speciesReference species="A" constant="true"/></listOfReactants>\n' ...
                   '<listOfModifiers><modifierSpeciesReference species="A"/></listOfModifiers>\n' ...
                   '<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>A</ci></math>\n' ...
                   '<listOfLocalParameters><localParameter id="k"/></listOfLocalParameters></kineticLaw>\n' ...
                   '</reaction></listOfReactions>\n' ...
                   '</model>\n</sbml>\n'], sbml_level, sbml_version, sbml_level, sbml_version);
  end

  filename = [tempname '.xml'];
  fid = fopen(filename,'w'); fprintf(fid,'%s',xml); fclose(fid);
  m = TranslateSBML(filename);
  delete(filename);

  t.compartment      = m.compartment;
  t.species          = m.species;
  t.parameter        = m.parameter;
  t.speciesReference = m.reaction.reactant;
  t.modifier         = m.reaction.modifier;
  if isfield(t.speciesReference,'isSetStoichiometry')  % SBML level 3
    t.speciesReference.isSetStoichiometry = 1;
  end

  t.kineticLaw = m.reaction.kineticLaw;
  if isfield(t.kineticLaw,'localParameter')
    t.kinetic_law_parameter_field = 'localParameter';
  else
    t.kinetic_law_parameter_field = 'parameter';
  end
  t.kinetic_law_parameter = t.kineticLaw.(t.kinetic_law_parameter_field);
  t.kineticLaw.(t.kinetic_law_parameter_field) = t.kinetic_law_parameter([]);
  t.kineticLaw.formula = '';
  t.kineticLaw.math    = '';

  t.reaction            = m.reaction;
  t.reaction.reactant   = t.speciesReference([]);
  t.reaction.product    = t.speciesReference([]);
  t.reaction.modifier   = t.modifier([]);
  t.reaction.kineticLaw = [];

  t.model             = m;
  t.model.id          = '';
  t.model.compartment = t.compartment([]);
  t.model.species     = t.species([]);
  t.model.parameter   = t.parameter([]);
  t.model.reaction    = t.reaction([]);


  % -------------------------------------------------------------------

function kinetic_law = make_kinetic_law(t, formula, par_ids, par_values)

  kinetic_law         = t.kineticLaw;
  kinetic_law.formula = formula;
  kinetic_law.math    = formula;
  for it = 1:length(par_ids)
    kinetic_law.(t.kinetic_law_parameter_field)(it) = make_parameter(t.kinetic_law_parameter, par_ids{it}, par_values(it));
  end


  % -------------------------------------------------------------------

function parameter = make_parameter(template, id, value)

  parameter            = template;
  parameter.id         = id;
  parameter.name       = id;
  parameter.value      = value;
  parameter.isSetValue = 1;
  if isfield(parameter,'constant')
    parameter.constant = 1;
  end


  % -------------------------------------------------------------------

function species = set_initial_concentration(species, value)

  species.initialConcentration      = value;
  species.isSetInitialConcentration = double(isfinite(value));
