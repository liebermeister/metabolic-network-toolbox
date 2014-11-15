function sbtab_document = network_to_sbtab(network,options)

% sbtab_document = network_to_sbtab(network,options)
%
% translate 'network' datastructure into 'sbtab' datastructure
% 
% If options.filename is given -> write output files: 
%  [filename]_Reaction.csv
%  [filename]_Compound.csv
%  [filename]_QuantityData.csv (only if argument modular_rate_law_table==1)
%
% Options are given in the structure 'options' with fields:
%
% only_reaction_table   : (default 0) produce only 'Reactions' table
% modular_rate_law_table: produce table with modular rate law kinetic constants 
% filename              : file name for saving the sbtab (extension is added automatically)
% save_in_one_file      : flag for saving SBtab in one file (default = 0)

try
  sbtab_version;
catch err,
  error('The SBtab toolbox for matlab must be installed');
end

eval(default('options','struct'));
options_default = struct('filename',[],'only_reaction_table',0,'modular_rate_law_table',1,'use_sbml_ids',0,'verbose',1,'write_concentrations',1,'save_in_one_file',0);
options         = join_struct(options_default,options);

if ~isfield(network,'kinetics'),
  options.modular_rate_law_table = 0;
end

% if necessary, make metabolite + reaction names syntactically correct:

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);
network.formulae = network_print_formulae(network,network.actions,network.metabolites);

reaction_table = sbtab_table_construct(struct('TableName','Reaction','TableType','Reaction'),{'Reaction','SumFormula'},{network.actions,network.formulae});

if isfield(network, 'reaction_names'),
  reaction_table = sbtab_table_add_column(reaction_table,'Name', network.reaction_names);
end

if isfield(network, 'reaction_KEGGID'),
  reaction_table = sbtab_table_add_column(reaction_table,'Reaction:Identifiers:kegg.reaction', network.reaction_KEGGID);
end

if isfield(network, 'reversible'),
  reaction_table = sbtab_table_add_column(reaction_table,'IsReversible',network.reversible);
end

if isfield(network, 'genes'),
  reaction_table = sbtab_table_add_column(reaction_table,'Gene',network.genes);
end

% metabolic regulation of enzymes
if norm(network.regulation_matrix)>0,
[nm,nr] = size(network.N);
for it = 1:nr,
  my_line = '';
  ind_pos = find(network.regulation_matrix(it,:)>0);
  ind_neg = find(network.regulation_matrix(it,:)<0);
  for itt = 1:length(ind_pos), 
    my_line = [my_line '+ ' network.metabolites{ind_pos(itt)} ' '];
  end
  for itt = 1:length(ind_neg), 
    my_line = [my_line '- ' network.metabolites{ind_neg(itt)} ' '];
  end
  if length(my_line) , my_line = my_line(1:end-1); end
  network.MetabolicRegulation{it,1} = my_line;
end
end

if isfield(network, 'TranscriptionalRegulation'),
  reaction_table = sbtab_table_add_column(reaction_table,'TranscriptionalRegulation',network.TranscriptionalRegulation);
end

if isfield(network, 'EC'),
  reaction_table = sbtab_table_add_column(reaction_table,'Enzyme:Identifiers:ec-code',network.EC);
end

if isfield(network, 'sbml_id_reaction'),
  reaction_table = sbtab_table_add_column(reaction_table,'SBML:reaction:id',network.sbml_id_reaction);
end

if isfield(network,'MetabolicRegulation'),
  reaction_table = sbtab_table_add_column(reaction_table,'MetabolicRegulation',network.MetabolicRegulation);
end

sbtab_document = sbtab_document_construct(struct,{'Reaction'},{reaction_table});;

if ~options.only_reaction_table,

  compound_table = sbtab_table_construct(struct('TableName','Compound','TableType','Compound'),{'Compound'},{network.metabolites});

  if isfield(network, 'metabolite_names'),
    compound_table = sbtab_table_add_column(compound_table,'Name', network.metabolite_names);
  end

  if isfield(network, 'metabolite_KEGGID'),
    compound_table = sbtab_table_add_column(compound_table,'Compound:Identifiers:kegg.compound',network.metabolite_KEGGID);
  end
  
  if isfield(network, 'sbml_id_species'),
    compound_table = sbtab_table_add_column(compound_table,'SBML:species:id',network.sbml_id_species);
  end

  compound_table = sbtab_table_add_column(compound_table,'External', uint8(network.external));

  if isfield(network, 'is_cofactor'),
    compound_table = sbtab_table_add_column(compound_table,'IsCofactor',network.is_cofactor);
  end

  sbtab_document = sbtab_document_add_table(sbtab_document,'Compound',compound_table);
end

if options.modular_rate_law_table,
  quantity_table = modular_rate_law_to_sbtab(network,[],struct('use_sbml_ids',options.use_sbml_ids,'write_concentrations',options.write_concentrations));
  sbtab_document = sbtab_document_add_table(sbtab_document,'Quantity',quantity_table);
end

if ~isempty(options.filename),
  switch options.save_in_one_file,
    case 0, sbtab_document_save(sbtab_document,options.filename,0,options.verbose);
    case 1, 
      if ~strcmp(options.filename(end-3:end),'.csv'),
        options.filename = [options.filename, '.csv'];
        end
      sbtab_document_save_to_one(sbtab_document,options.filename );
  end
end
