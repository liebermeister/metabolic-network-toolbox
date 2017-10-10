function quantity_table = modular_rate_law_to_sbtab(network,filename,options)

% quantity_table = modular_rate_law_to_sbtab(network,filename,options)
% 
% Options and default values:  
%
%   options.use_sbml_ids = 1;
%   options.write_individual_kcat = 1;
%   options.write_std = 1;
%   options.write_concentrations = 1;
%   options.write_enzyme_concentrations = 1;
%   options.kinetics_mean = [];
%   options.kinetics_std = [];
%   options.modular_rate_law_parameter_id = 0;
%   options.document_name = '';
%
% Requires the SBtab toolbox

eval(default('options','struct','filename','[]'));

options_default = struct('use_sbml_ids',1,'write_individual_kcat',1,'write_std',1, 'write_concentrations',1,'write_enzyme_concentrations',1,'modular_rate_law_parameter_id',0,'document_name','','kinetics_mean', [], 'kinetics_std', []);
options         = join_struct(options_default,options);

kinetics_mean = options.kinetics_mean;
kinetics_std  = options.kinetics_std;

if ~isfield(network.kinetics,'u'),
  options.write_enzyme_concentrations = 0;
elseif isempty(network.kinetics.u),
  options.write_enzyme_concentrations = 0;
end

if ~isfield(network.kinetics,'c'),
  options.write_concentrations = 0;
elseif isempty(network.kinetics.c),
  options.write_concentrations = 0;
end

switch network.kinetics.type
  case {'cs','ds','ms','rp','ma','fm'}, % UPDATE rate law names!
  otherwise, error('Conversion is only possible for modular rate law');
end

[nm,nr] = size(network.N);

if options.use_sbml_ids, 
  metabolites = network.sbml_id_species; 
  reactions   = network.sbml_id_reaction; 
else
  metabolites = network.metabolites; 
  reactions   = network.actions; 
  [metabolites,reactions] = network_adjust_names_for_sbml_export(metabolites,reactions);
end

if isfield(network, 'metabolite_KEGGID'),
  metabolite_KEGGID = network.metabolite_KEGGID; 
else
  metabolite_KEGGID = repmat({''},nm,1);
end

if isfield(network, 'reaction_KEGGID'),
  reaction_KEGGID = network.reaction_KEGGID; 
else
  reaction_KEGGID = repmat({''},nr,1);
end

[nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(network);

column_quantity = [ repmat({'equilibrium constant'},nr,1); ...
                    repmat({'catalytic rate constant geometric mean'},nr,1); ...
                    repmat({'Michaelis constant'},nKM,1); ...
                    repmat({'activation constant'},nKA,1); ...
                    repmat({'inhibitory constant'},nKI,1); ...
                   ];

if ~isfield(network.kinetics,'KV'),
  network.kinetics.KV = sqrt(network.kinetics.Kcatf .* network.kinetics.Kcatr);
  if length(kinetics_mean), 
    %% thius could also be taken from the actual results ..
    kinetics_mean.KV = nan * network.kinetics.KV;
    kinetics_std.KV  = nan * network.kinetics.KV;
  end
end

column_median = full( [ network.kinetics.Keq; ...
                    network.kinetics.KV; ...
                    column(network.kinetics.KM(KM_indices)); ...
                    column(network.kinetics.KA(KA_indices)); ...
                    column(network.kinetics.KI(KI_indices)); ...
                   ]);

if length(kinetics_mean), 
  column_mean = full( [kinetics_mean.Keq; ...
                      kinetics_mean.KV; ...
                      column(kinetics_mean.KM(KM_indices)); ...
                      column(kinetics_mean.KA(KA_indices)); ...
                      column(kinetics_mean.KI(KI_indices)); ...
                      ]);

  column_std = full( [kinetics_std.Keq; ...
                      kinetics_std.KV; ...
                      column(kinetics_std.KM(KM_indices)); ...
                      column(kinetics_std.KA(KA_indices)); ...
                      column(kinetics_std.KI(KI_indices)); ...
                      ]);
end

column_unit = [ repmat({'dimensionless'},nr,1); ...
                repmat({'1/s'},nr,1); ...
                repmat({'mM'},nKM,1); ...
                repmat({'mM'},nKA,1); ...
                repmat({'mM'},nKI,1); ...
              ];    

[iKM,jKM] = ind2sub([nr,nm],KM_indices);
[iKA,jKA] = ind2sub([nr,nm],KA_indices);
[iKI,jKI] = ind2sub([nr,nm],KI_indices);


dumKM = {};
for it = 1:length(iKM),
  dumKM{it,1} = sprintf('kM_R%d_%s', iKM(it), network.metabolites{jKM(it)});
end

dumKA = {};
for it = 1:length(iKA),
  dumKA{it,1} = sprintf('kA_R%d_%s', iKA(it), network.metabolites{jKA(it)});
end

dumKI = {};
for it = 1:length(iKI),
  dumKI{it,1} = sprintf('kI_R%d_%s', iKI(it), network.metabolites{jKI(it)});
end


column_parameterID = [numbered_names_simple('kEQ_R',nr); ...
                    numbered_names_simple('kC_R',nr); ...
                    dumKM; ...
                    dumKA; ...
                    dumKI; ...
                   ];


column_reaction = [ reactions; ...
                    reactions; ...
                    reactions(iKM); ...
                    reactions(iKA); ...
                    reactions(iKI); ...
                   ];    

column_reaction_KEGGID = [ reaction_KEGGID; ...
                    reaction_KEGGID; ...
                    reaction_KEGGID(iKM); ...
                    reaction_KEGGID(iKA); ...
                    reaction_KEGGID(iKI); ...
                   ];    

column_compound = [ repmat({''},nr,1); ...
                    repmat({''},nr,1); ...
                    metabolites(jKM); ...
                    metabolites(jKA); ...
                    metabolites(jKI); ...
                  ];

column_compound_KEGGID = [ repmat({''},nr,1); ...
                    repmat({''},nr,1); ...
                    metabolite_KEGGID(jKM); ...
                    metabolite_KEGGID(jKA); ...
                    metabolite_KEGGID(jKI); ...
                  ];

% Only if required: write individual kcat values

if ~isfield(network.kinetics,'Kcatf'),
  options.write_individual_kcat = 0;
end

if options.write_individual_kcat,

column_quantity = [ column_quantity; ...
                    repmat({'substrate catalytic rate constant'},nr,1); ...
                    repmat({'product catalytic rate constant'},nr,1); ...
                  ];

column_parameterID = [ column_parameterID; ...
                    numbered_names_simple('kcrf_R',nr); ...
                    numbered_names_simple('kcrr_R',nr); ...
                   ];

column_median = [ column_median; ...
                 network.kinetics.Kcatf; ...
                 network.kinetics.Kcatr; ...
               ];

if length(kinetics_mean), 
column_mean = [ column_mean; ...
                 kinetics_mean.Kcatf; ...
                 kinetics_mean.Kcatr; ...
               ];

column_std = [ column_std; ...
                 kinetics_std.Kcatf; ...
                 kinetics_std.Kcatr; ...
               ];
end

column_unit = [ column_unit; ...
                repmat({'1/s'},nr,1); ...
                repmat({'1/s'},nr,1); ...
              ];    

column_reaction = [ column_reaction; ...
                    reactions; ...
                    reactions; ...
                   ];    

column_compound = [ column_compound; ...
                    repmat({''},nr,1); ...
                    repmat({''},nr,1); ...
                  ];

column_reaction_KEGGID = [ column_reaction_KEGGID; ...
                    reaction_KEGGID; ...
                    reaction_KEGGID; ...
                   ];    

column_compound_KEGGID = [ column_compound_KEGGID; ...
                    repmat({''},nr,1); ...
                    repmat({''},nr,1); ...
                  ];
  
end


% Only if required: write metabolite and enzyme concentrations


if options.write_concentrations,
  column_quantity        = [ column_quantity;  repmat({'concentration'},nm,1)];
  column_parameterID     = [ column_parameterID; cellstr([char(repmat({'c_'},nm,1)) char(network.metabolites)])];
  column_median          = [ column_median; network.kinetics.c];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; kinetics_mean.c];
    column_std           = [ column_std; kinetics_std.c];
  end
  column_unit            = [ column_unit;  repmat({'mM'},nm,1)];    
  column_reaction        = [ column_reaction;  repmat({''},nm,1)];    
  column_compound        = [ column_compound; metabolites];
  column_reaction_KEGGID = [ column_reaction_KEGGID;  repmat({''},nm,1)];    
  column_compound_KEGGID = [ column_compound_KEGGID;  metabolite_KEGGID];
end


if options.write_enzyme_concentrations,
  column_quantity        = [ column_quantity; repmat({'concentration of enzyme'},nr,1)];
  column_parameterID     = [ column_parameterID; numbered_names_simple('u_R',nr)];
  column_median           = [ column_median; network.kinetics.u];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; kinetics_mean.u];
    column_std           = [ column_std; kinetics_std.u];
  end
  column_unit            = [ column_unit; repmat({'mM'},nr,1)];
  column_reaction        = [ column_reaction;  reactions];    
  column_compound        = [ column_compound; repmat({''},nr,1)];
  column_reaction_KEGGID = [ column_reaction_KEGGID; reaction_KEGGID];    
  column_compound_KEGGID = [ column_compound_KEGGID;  repmat({''},nr,1)];
end


if isfield(network,'metabolite_mass'),
  column_quantity        = [ column_quantity;  repmat({'molecular mass'},nm,1)];
  column_parameterID     = [ column_parameterID; cellstr([char(repmat({'mm_'},nm,1)) char(network.metabolites)])];
  column_median           = [ column_median; network.metabolite_mass];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; repmat({''},nm,1)];
    column_std           = [ column_std; repmat({''},nm,1)];
  end
  column_unit            = [ column_unit;  repmat({'Da'},nm,1)];    
  column_reaction        = [ column_reaction;  repmat({''},nm,1)];
  column_compound        = [ column_compound; metabolites];
  column_reaction_KEGGID = [ column_reaction_KEGGID;  repmat({''},nm,1)];    
  column_compound_KEGGID = [ column_compound_KEGGID;  metabolite_KEGGID];
end

if isfield(network,'enzyme_mass'),
  column_quantity        = [ column_quantity; repmat({'protein molecular mass'},nr,1)];
  column_parameterID     = [ column_parameterID; numbered_names_simple('em_',nr)];
  column_median           = [ column_median; network.enzyme_mass];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; repmat({''},nr,1)];
    column_std           = [ column_std; repmat({''},nr,1)];
  end
  column_unit            = [ column_unit; repmat({'Da'},nr,1)];
  column_reaction        = [ column_reaction;  reactions];    
  column_compound        = [ column_compound; repmat({''},nr,1)];
  column_reaction_KEGGID = [ column_reaction_KEGGID; reaction_KEGGID];    
  column_compound_KEGGID = [ column_compound_KEGGID;  repmat({''},nr,1)];
end

if isfield(network, 'metabolite_KEGGID'),
  if isfield(network, 'reaction_KEGGID'), 
      if length(kinetics_mean), 
        quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit','Mean','Std','Reaction:Identifiers:kegg.reaction','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_median,column_unit,column_mean,column_std,column_reaction_KEGGID,column_compound_KEGGID});
      else,
        quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit','Reaction:Identifiers:kegg.reaction','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_median,column_unit,column_reaction_KEGGID,column_compound_KEGGID});
      end
      else
        if length(kinetics_mean), 
        quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit','Mean','Std','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_median,column_unit,column_mean,column_std,column_compound_KEGGID});
      else,
        quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_median,column_unit,column_compound_KEGGID});
      end
  end
else,
  if isfield(network, 'reaction_KEGGID'),  
        if length(kinetics_mean), 
    quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit','Mean','Std','Reaction:Identifiers:kegg.reaction'}, {column_quantity,column_reaction,column_compound,column_median,column_unit,column_mean,column_std,column_reaction_KEGGID});        
        else  
    quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit','Reaction:Identifiers:kegg.reaction'}, {column_quantity,column_reaction,column_compound,column_median,column_unit,column_reaction_KEGGID});        
          end
  else
        if length(kinetics_mean), 
    quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit','Mean','Std'}, {column_quantity,column_reaction,column_compound,column_median,column_unit,column_mean,column_std});
        else
    quantity_table = sbtab_table_construct(struct('TableName','Parameter','TableType','Quantity','Document',options.document_name), {'QuantityType','Reaction','Compound','Median','Unit'}, {column_quantity,column_reaction,column_compound,column_median,column_unit});
        end
        end
end


if options.modular_rate_law_parameter_id, 
  quantity_table = sbtab_table_add_column(quantity_table, 'ID', column_parameterID, 1);
end

if length(filename), 
  sbtab_table_save(quantity_table,struct('filename',filename)); 
end
