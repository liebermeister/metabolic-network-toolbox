function quantity_table = modular_rate_law_to_sbtab(network,filename,options)

% quantity_table = modular_rate_law_to_sbtab(network,filename,options)
% 
% Options and default values:  
%
%   options.use_sbml_ids = 1;
%   options.write_std = 1;
%   options.write_all_quantities        = 'no';
%   options.write_individual_kcat       = 1;
%   options.write_concentrations        = 1;
%   options.write_enzyme_concentrations = 1;
%   options.write_standard_chemical_potential = 0;
%   options.write_chemical_potential    = 0;
%   options.write_reaction_affinity     = 0;
%   options.kinetics_mean = [];
%   options.kinetics_mode = [];
%   options.kinetics_std = [];
%   options.kinetics_geom_std = [];
%   options.kinetics_geom_mean = [];
%   options.modular_rate_law_parameter_id = 0;
%   options.document_name = '';
%
% Requires the SBtab toolbox

eval(default('options','struct','filename','[]'));

options_default = struct('use_sbml_ids',1,'write_individual_kcat',1,'write_std',1, 'write_concentrations',1,'write_enzyme_concentrations',1,'write_all_quantities', 'no', 'write_standard_chemical_potential', 0, 'write_chemical_potential', 0, 'write_std_GFE_of_reaction', 0, 'write_reaction_affinity', 0, 'modular_rate_law_parameter_id',0,'document_name','','kinetics_mean', [], 'kinetics_std', [],'kinetics_mode', [],  'kinetics_geom_mean', [],'kinetics_geom_std', []);
options         = join_struct(options_default,options);

switch options.write_all_quantities,
  case 'all',
   options.write_individual_kcat       = 1;
   options.write_concentrations        = 1;
   options.write_enzyme_concentrations = 1;
   options.write_standard_chemical_potential = 1;
   options.write_chemical_potential    = 1;
   options.write_reaction_affinity     = 1;
   options.write_std_GFE_of_reaction   = 1;
  case 'many',
   options.write_individual_kcat       = 1;
   options.write_concentrations        = 1;
   options.write_enzyme_concentrations = 1;
   options.write_standard_chemical_potential = 1;
   options.write_chemical_potential    = 0;
   options.write_reaction_affinity     = 0;
   options.write_std_GFE_of_reaction   = 0;
end

kinetics_mode = options.kinetics_mode;
kinetics_mean = options.kinetics_mean;
kinetics_std  = options.kinetics_std;
kinetics_geom_mean  = options.kinetics_geom_mean;
kinetics_geom_std  = options.kinetics_geom_std;

if ~isfield(network.kinetics,'Kcatf'),
  options.write_individual_kcat = 0;
end

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

if ~isfield(network.kinetics,'mu0'),
  options.write_standard_chemical_potential = 0;
elseif isempty(network.kinetics.mu0),
  options.write_standard_chemical_potential = 0;
end

if ~isfield(network.kinetics,'mu'),
  options.write_chemical_potential = 0;
elseif isempty(network.kinetics.mu),
  options.write_chemical_potential = 0;
end

if ~isfield(network.kinetics,'A'),
  options.write_reaction_affinity = 0;
elseif isempty(network.kinetics.A),
  options.write_reaction_affinity = 0;
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
    %% this could also be taken from the actual results ..
    kinetics_mean.KV = nan * network.kinetics.KV;
    kinetics_std.KV  = nan * network.kinetics.KV;
    kinetics_mode.KV = nan * network.kinetics.KV;
    kinetics_geom_mean.KV  = nan * network.kinetics.KV;
    kinetics_geom_std.KV  = nan * network.kinetics.KV;
  end
end

column_mode = full([ network.kinetics.Keq; ...
                    network.kinetics.KV; ...
                    column(network.kinetics.KM(KM_indices)); ...
                    column(network.kinetics.KA(KA_indices)); ...
                    column(network.kinetics.KI(KI_indices)); ...
                   ]);

column_unit = [ repmat({'dimensionless'},nr,1); ...
                repmat({'1/s'},nr,1); ...
                repmat({'mM'},nKM,1); ...
                repmat({'mM'},nKA,1); ...
                repmat({'mM'},nKI,1); ...
              ];    

if length(kinetics_mean), 
  column_mean = full( [kinetics_mean.Keq; ...
                      kinetics_mean.KV; ...
                      column(kinetics_mean.KM(KM_indices)); ...
                      column(kinetics_mean.KA(KA_indices)); ...
                      column(kinetics_mean.KI(KI_indices)); ...
                      ]);
end

if length(kinetics_std), 
  column_std = full( [kinetics_std.Keq; ...
                      kinetics_std.KV; ...
                      column(kinetics_std.KM(KM_indices)); ...
                      column(kinetics_std.KA(KA_indices)); ...
                      column(kinetics_std.KI(KI_indices)); ...
                      ]);
end

if length(kinetics_geom_mean), 
  column_geom_mean = full( [kinetics_geom_mean.Keq; ...
                      kinetics_geom_mean.KV; ...
                      column(kinetics_geom_mean.KM(KM_indices)); ...
                      column(kinetics_geom_mean.KA(KA_indices)); ...
                      column(kinetics_geom_mean.KI(KI_indices)); ...
                      ]);
end

if length(kinetics_geom_std), 
  column_geom_std = full( [kinetics_geom_std.Keq; ...
                      kinetics_geom_std.KV; ...
                      column(kinetics_geom_std.KM(KM_indices)); ...
                      column(kinetics_geom_std.KA(KA_indices)); ...
                      column(kinetics_geom_std.KI(KI_indices)); ...
                      ]);
end


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

if options.write_individual_kcat,

column_quantity = [ column_quantity; ...
                    repmat({'substrate catalytic rate constant'},nr,1); ...
                    repmat({'product catalytic rate constant'},nr,1); ...
                  ];

column_parameterID = [ column_parameterID; ...
                    numbered_names_simple('kcrf_R',nr); ...
                    numbered_names_simple('kcrr_R',nr); ...
                   ];

column_mode = [ column_mode; ...
                 network.kinetics.Kcatf; ...
                 network.kinetics.Kcatr; ...
               ];

if length(kinetics_mean), 
column_mean = [ column_mean; ...
                 kinetics_mean.Kcatf; ...
                 kinetics_mean.Kcatr; ...
               ];
end

if length(kinetics_std), 
column_std = [ column_std; ...
                 kinetics_std.Kcatf; ...
                 kinetics_std.Kcatr; ...
               ];
end

if length(kinetics_geom_mean), 
column_geom_mean = [ column_geom_mean; ...
                 kinetics_geom_mean.Kcatf; ...
                 kinetics_geom_mean.Kcatr; ...
               ];
end

if length(kinetics_geom_std), 
column_geom_std = [ column_geom_std; ...
                 kinetics_geom_std.Kcatf; ...
                 kinetics_geom_std.Kcatr; ...
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


if options.write_standard_chemical_potential,
  column_quantity        = [ column_quantity; repmat({'standard chemical potential'},nm,1)];
  column_parameterID     = [ column_parameterID; numbered_names_simple('mu0_M',nm)];
  column_mode           = [ column_mode; network.kinetics.mu0];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; kinetics_mean.mu0];
    column_std           = [ column_std; kinetics_std.mu0];
    column_geom_mean           = [ column_geom_mean; kinetics_geom_mean.mu0];
    column_geom_std           = [ column_geom_std; kinetics_geom_std.mu0];
  end
  column_unit            = [ column_unit; repmat({'kJ/mol'},nm,1)];
  column_reaction        = [ column_reaction; repmat({''},nm,1)];    
  column_compound        = [ column_compound; metabolites];
  column_reaction_KEGGID = [ column_reaction_KEGGID;  repmat({''},nm,1)];    
  column_compound_KEGGID = [ column_compound_KEGGID; metabolite_KEGGID];
end

if options.write_chemical_potential,
  column_quantity        = [ column_quantity; repmat({'chemical_potential'},nm,1)];
  column_parameterID     = [ column_parameterID; numbered_names_simple('mu_M',nm)];
  column_mode           = [ column_mode; network.kinetics.mu];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; kinetics_mean.mu];
    column_std           = [ column_std; kinetics_std.mu];
    column_geom_mean           = [ column_geom_mean; kinetics_geom_mean.mu];
    column_geom_std           = [ column_geom_std; kinetics_geom_std.mu];
  end
  column_unit            = [ column_unit; repmat({'kJ/mol'},nm,1)];
  column_reaction        = [ column_reaction; repmat({''},nm,1)];    
  column_compound        = [ column_compound; metabolites];
  column_reaction_KEGGID = [ column_reaction_KEGGID;  repmat({''},nm,1)];    
  column_compound_KEGGID = [ column_compound_KEGGID; metabolite_KEGGID];
end

if options.write_reaction_affinity,
  column_quantity        = [ column_quantity; repmat({'reaction affinity'},nr,1)];
  column_parameterID     = [ column_parameterID; numbered_names_simple('A_R',nr)];
  column_mode           = [ column_mode; network.kinetics.A];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; kinetics_mean.A];
    column_std           = [ column_std; kinetics_std.A];
    column_geom_mean      = [ column_geom_mean; kinetics_geom_mean.A];
    column_geom_std      = [ column_geom_std; kinetics_geom_std.A];
  end
  column_unit            = [ column_unit; repmat({'kJ/mol'},nr,1)];
  column_reaction        = [ column_reaction;  reactions];    
  column_compound        = [ column_compound; repmat({''},nr,1)];
  column_reaction_KEGGID = [ column_reaction_KEGGID; reaction_KEGGID];    
  column_compound_KEGGID = [ column_compound_KEGGID;  repmat({''},nr,1)];
end

if options.write_std_GFE_of_reaction,
  column_quantity        = [ column_quantity; repmat({'standard Gibbs energy of reaction'},nr,1)];
  column_parameterID     = [ column_parameterID; numbered_names_simple('deltaGFE0_R',nr)];
  column_mode            = [ column_mode; -RT * log(network.kinetics.Keq)];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; -RT * log(kinetics_mean.Keq)];
    column_std           = [ column_std; nan * kinetics_std.A];
    column_geom_mean     = [ column_geom_mean; nan * kinetics_geom_mean.A];
    column_geom_std      = [ column_geom_std; nan * kinetics_geom_std.A];
  end
  column_unit            = [ column_unit; repmat({'kJ/mol'},nr,1)];
  column_reaction        = [ column_reaction;  reactions];    
  column_compound        = [ column_compound; repmat({''},nr,1)];
  column_reaction_KEGGID = [ column_reaction_KEGGID; reaction_KEGGID];    
  column_compound_KEGGID = [ column_compound_KEGGID;  repmat({''},nr,1)];
end

if options.write_concentrations,
  column_quantity        = [ column_quantity;  repmat({'concentration'},nm,1)];
  column_parameterID     = [ column_parameterID; cellstr([char(repmat({'c_'},nm,1)) char(network.metabolites)])];
  column_mode          = [ column_mode; network.kinetics.c];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; kinetics_mean.c];
    column_std           = [ column_std; kinetics_std.c];
    column_geom_mean           = [ column_geom_mean; kinetics_geom_mean.c];
    column_geom_std           = [ column_geom_std; kinetics_geom_std.c];
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
  column_mode           = [ column_mode; network.kinetics.u];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; kinetics_mean.u];
    column_std           = [ column_std; kinetics_std.u];
    column_geom_mean           = [ column_geom_mean; kinetics_geom_mean.u];
    column_geom_std           = [ column_geom_std; kinetics_geom_std.u];
  end
  column_unit            = [ column_unit; repmat({'mM'},nr,1)];
  column_reaction        = [ column_reaction;  reactions];    
  column_compound        = [ column_compound; repmat({''},nr,1)];
  column_reaction_KEGGID = [ column_reaction_KEGGID; reaction_KEGGID];    
  column_compound_KEGGID = [ column_compound_KEGGID;  repmat({''},nr,1)];
end

if isfield(network,'enzyme_mass'),
  column_quantity        = [ column_quantity; repmat({'protein molecular mass'},nr,1)];
  column_parameterID     = [ column_parameterID; numbered_names_simple('em_',nr)];
  column_mode           = [ column_mode; network.enzyme_mass];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; repmat({''},nr,1)];
    column_std           = [ column_std; repmat({''},nr,1)];
    column_geom_mean           = [ column_geom_mean; repmat({''},nr,1)];
    column_geom_std           = [ column_geom_std; repmat({''},nr,1)];
  end
  column_unit            = [ column_unit; repmat({'Da'},nr,1)];
  column_reaction        = [ column_reaction;  reactions];    
  column_compound        = [ column_compound; repmat({''},nr,1)];
  column_reaction_KEGGID = [ column_reaction_KEGGID; reaction_KEGGID];    
  column_compound_KEGGID = [ column_compound_KEGGID;  repmat({''},nr,1)];
end

if isfield(network,'metabolite_mass'),
  column_quantity        = [ column_quantity;  repmat({'molecular mass'},nm,1)];
  column_parameterID     = [ column_parameterID; cellstr([char(repmat({'mm_'},nm,1)) char(network.metabolites)])];
  column_mode           = [ column_mode; network.metabolite_mass];
  if length(kinetics_mean), 
    column_mean          = [ column_mean; repmat({''},nm,1)];
    column_std           = [ column_std; repmat({''},nm,1)];
    column_geom_mean      = [ column_geom_mean; repmat({''},nm,1)];
    column_geom_std      = [ column_geom_std; repmat({''},nm,1)];
  end
  column_unit            = [ column_unit;  repmat({'Da'},nm,1)];    
  column_reaction        = [ column_reaction;  repmat({''},nm,1)];
  column_compound        = [ column_compound; metabolites];
  column_reaction_KEGGID = [ column_reaction_KEGGID;  repmat({''},nm,1)];    
  column_compound_KEGGID = [ column_compound_KEGGID;  metabolite_KEGGID];
end


sbtab_attributes = struct('SBtabVersion',num2str(sbtab_version),'TableType','Quantity','TableName','Parameter','Document',options.document_name);

if isfield(network, 'metabolite_KEGGID'),
  if isfield(network, 'reaction_KEGGID'), 
      if length(kinetics_mean), 
        quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Mode','Unit','UnconstrainedGeometricMean','UnconstrainedGeometricStd','UnconstrainedMean','UnconstrainedStd','Reaction:Identifiers:kegg.reaction','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_mode,column_unit,column_geom_mean,column_geom_std,column_mean,column_std,column_reaction_KEGGID,column_compound_KEGGID});
      else,
        quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Value','Unit','Reaction:Identifiers:kegg.reaction','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_mode,column_unit,column_reaction_KEGGID,column_compound_KEGGID});
      end
      else
        if length(kinetics_mean), 
        quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Mode','Unit','UnconstrainedGeometricMean','UnconstrainedGeometricStd','UnconstrainedMean','UnconstrainedStd','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_mode,column_unit,column_geom_mean,column_geom_std,column_mean,column_std,column_compound_KEGGID});
      else,
        quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Value','Unit','Compound:Identifiers:kegg.compound'}, {column_quantity,column_reaction,column_compound,column_mode,column_unit,column_compound_KEGGID});
      end
  end
else,
  if isfield(network, 'reaction_KEGGID'),  
        if length(kinetics_mean), 
    quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Mode','Unit','UnconstrainedGeometricMean','UnconstrainedGeometricStd','UnconstrainedMean','UnconstrainedStd','Reaction:Identifiers:kegg.reaction'}, {column_quantity,column_reaction,column_compound,column_mode,column_unit,column_geom_mean,column_geom_std,column_mean,column_std,column_reaction_KEGGID});        
        else  
    quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Value','Unit','Reaction:Identifiers:kegg.reaction'}, {column_quantity,column_reaction,column_compound,column_mode,column_unit,column_reaction_KEGGID});        
          end
  else
        if length(kinetics_mean), 
    quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Mode','Unit','UnconstrainedGeometricMean','UnconstrainedGeometricStd'},'UnconstrainedMean','UnconstrainedStd', {column_quantity,column_reaction,column_compound,column_mode,column_unit,column_geom_mean,column_geom_std,column_mean,column_std});
        else
    quantity_table = sbtab_table_construct(sbtab_attributes, {'QuantityType','Reaction','Compound','Value','Unit'}, {column_quantity,column_reaction,column_compound,column_mode,column_unit});
        end
        end
end


if options.modular_rate_law_parameter_id, 
  quantity_table = sbtab_table_add_column(quantity_table, 'ID', column_parameterID, 1);
end

if length(filename), 
  sbtab_table_save(quantity_table,struct('filename',filename)); 
end
