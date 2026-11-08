function [network, kinetic_data] = sbtab_to_network(filename,options)

% SBTAB_TO_NETWORK Read a network structure and kinetic data from SBtab
%
% [network, parameter_table] = sbtab_to_network(filename,options)
%
% Read network structure and kinetic data from SBtab files
% (for the opposite conversion, use 'network_to_sbtab');
%
% The second output contains the SBtab table with kinetic parameters
% (empty if no parameters were read).
%
%Insert kinetic parameters that *directly* appear in the rate law specified in "options.kinetic_law"
%
% Assumes three input files:
%  [filename]_Reaction.tsv
%  [filename]_Compound.tsv
%  [filename]_QuantityData.tsv
%
% The last one can be omitted if options.load_quantity_table = 0;
% or another file name can be explicitly given in options.quantity_table_file
% also see network_to_sbtab.m
%
% Further optional input tables:
%  [filename]_Regulation.tsv  (table type Regulation; in an SBtab document: TableID 'Regulation')
%  [filename]_Regulator.tsv   (table type Regulator, regulators of gene expression with column
%                              ActivityCondition; in an SBtab document: TableID 'Regulator';
%                              see sbtab_regulator_table_to_network)
%  [filename]_Position.tsv (or [filename]_Layout.tsv)
%
% OR: instead of a filename, an SBtab data structure can given directly!!
%
% OR: if options.one_sbtab_file = 1 has been set, all data are read from one SBtab file (several tables)
%
% Options and their default values:
%
%  options_default.verbose             = 0;     % verbose mode
%  options_default.kinetic_law         = 'cs';  % type of modular rate law to be used
%  options_default.use_compound_table  = 1;     % use input file for compounds?
%  options_default.load_quantity_table = 1;     % use input file for parameter data?
%  options_default.quantity_table_file = '';    % filename for parameter data
%  options_default.use_sbml_ids        = 0;     % use sbml_ids from input files
%  options_default.extension           ='.tsv'; % default extension of input files
%  options_default.read_positions      = 1;     % expect SBtab table to contain position information
%  options_default.position_file       = [];    % filename for position information (overrides the default file name)
%  options_default.only_reaction_table = 0;     % read only reaction table (no compounds and parameters)
%                                               % assume full filename to be given
%  options_default.one_sbtab_file      = 1;     % all data from one SBtab file (several tables)
%
% See also: network_to_sbtab, network_sbml_import

try
  sbtab_version;
catch err
  error('The SBtab toolbox for matlab must be installed');
end

%network_sbtab_check_for_errors(my_sbtab);

eval(default('options','struct'));

options_default.verbose             = 0;     % verbose mode
options_default.kinetic_law         = 'cs';  % type of modular rate law to be used
options_default.use_compound_table  = 1;     % use input file for compounds?
options_default.load_quantity_table = 1;     % use input file for parameter data?
options_default.quantity_table_file = '';    % filename for parameter data
options_default.use_sbml_ids        = 0;     % use sbml_ids from input files
options_default.extension           ='.tsv'; % default extension of input files
options_default.read_positions      = 1;     % expect SBtab table to contain position information
options_default.position_file       = [];    % filename for position information
options_default.omit_water          = 1;     % set row for water in stoich matrix to zero
options_default.only_reaction_table = 0;     % read only reaction table (no compounds and parameters)
% assume full filename to be given
options_default.one_sbtab_file      = 1;     % all data from one SBtab file (several tables)

options = join_struct(options_default,options);

kinetic_data = [];

if options.only_reaction_table
  options.use_compound_table = 0;
  options.load_quantity_table = 0;
end

% sbtab_table_print(filename.tables.Parameter,'all')

if options.one_sbtab_file
  filename = sbtab_document_load_from_one(filename);
  if ~isfield(filename.tables,'Parameter'), options.load_quantity_table = 0; end
end

if ~length(options.position_file)
  if ~ischar(filename)
    if ~isfield(filename.tables,'Position')
      options.read_positions = 0;
    end
  end
end

if ischar(filename)
  if options.only_reaction_table
    filename_reactions = filename;
    filename_compounds = [];
    filename_regulation = [];
    filename_regulator = [];
    filename_position  = [];
  else
    % Assume that "filename" is an actual filename
    filename_reactions = [filename '_Reaction' options.extension];
    filename_compounds = [filename '_Compound' options.extension];
    filename_regulation = [filename '_Regulation' options.extension];
    if exist(filename_regulation,'file')
      filename_regulation = sbtab_table_load(filename_regulation);
    else
      filename_regulation = [];
    end
    filename_regulator = [filename '_Regulator' options.extension];
    if ~exist(filename_regulator,'file')
      filename_regulator = [];
    end
    if options.read_positions
      % layout table: file [filename]_Position or (newer name) [filename]_Layout
      filename_position = [filename '_Position' options.extension];
      if ~exist(filename_position,'file')
        filename_position = [filename '_Layout' options.extension];
      end
      if ~exist(filename_position,'file')
        filename_position = [];
      end
    else
      filename_position = [];
    end
  end
else
  % Assume that "filename" is an SBtab data structure!
  filename_reactions  = filename.tables.Reaction;
  if isfield(filename.tables,'Compound')
    filename_compounds  = filename.tables.Compound;
  else
    filename_compounds = [];
  end
  % Regulation table (kinetic regulation of enzymes): TableID 'Regulation'
  % (older files: 'EnzymeRegulation'), or another table of TableType 'Regulation'
  filename_regulation = [];
  if isfield(filename.tables,'Regulation')
    filename_regulation = filename.tables.Regulation;
  elseif isfield(filename.tables,'EnzymeRegulation')
    filename_regulation = filename.tables.EnzymeRegulation;
  else
    table_ids = fieldnames(filename.tables);
    for it = 1:length(table_ids)
      my_table = filename.tables.(table_ids{it});
      if isfield(my_table,'attributes') && isfield(my_table.attributes,'TableType') && strcmp(my_table.attributes.TableType,'Regulation')
        filename_regulation = my_table;
        break
      end
    end
  end
  % Regulator table (regulators of gene expression): TableID 'Regulator',
  % or another table of TableType 'Regulator'
  filename_regulator = [];
  if isfield(filename.tables,'Regulator')
    filename_regulator = filename.tables.Regulator;
  else
    table_ids = fieldnames(filename.tables);
    for it = 1:length(table_ids)
      my_table = filename.tables.(table_ids{it});
      if isfield(my_table,'attributes') && isfield(my_table.attributes,'TableType') && strcmp(my_table.attributes.TableType,'Regulator')
        filename_regulator = my_table;
        break
      end
    end
  end
  if isfield(filename.tables,'Position')
    filename_position = filename.tables.Position;
  elseif isfield(filename.tables,'Layout')
    filename_position = filename.tables.Layout;
  else
    filename_position = [];
  end
end

% a layout file given in options.position_file has priority
if length(options.position_file)
  filename_position = options.position_file;
end

% If thermodynamic parameters are in an extra table, copy them into table "Parameters"

if isstruct(filename) && isfield(filename.tables,'Thermodynamics')
  if isfield(filename.tables,'Parameter')
    filename.tables.Parameter = sbtab_table_combine_vertically(filename.tables.Thermodynamics,filename.tables.Parameter);
  end
end

if options.use_compound_table ==0
  filename_compounds = [];
end

network = network_build_from_sum_formulae(filename_reactions, filename_compounds, [], filename_regulation);

if length(filename_regulator)
  network = sbtab_regulator_table_to_network(network, filename_regulator);
end

if isfield(network,'IsReversible')
  network.reversible(find(strcmpi('True',network.IsReversible))) = 1;
  network.reversible(find(strcmpi('False',network.IsReversible))) = 0;
end

if isfield(network,'IsConstant')
  network.external(find(strcmpi('True',network.IsConstant))) = 1;
end

if options.omit_water
  if network_find_water(network)
    if options.verbose
      display('  WARNING (sbtab_to_network.m): Water as a chemical species causes problems in thermodynamic and kinetic models; setting stoichiometric coefficients for water to zero');
    end
    network.N(network_find_water(network),:) = 0;
  end
end


% -------------------------------------------
% Layout (optional)

if length(filename_position)
  network = netgraph_read_positions(network, filename_position);
end

% -------------------------------------------
% Parameters (optional)

if options.load_quantity_table

  if length(options.quantity_table_file)
    quantity_table = options.quantity_table_file;
  elseif ischar(filename)
    quantity_table = [filename '_QuantityData' options.extension];
  else
    quantity_table = filename.tables.Parameter;
  end
  [network.kinetics, kinetic_data, other_parameters] = sbtab_to_modular_rate_law(network,quantity_table,options);

  if sum(isfinite(other_parameters.metabolite_mass))
    network.metabolite_mass = other_parameters.metabolite_mass;
  end

  if sum(isfinite(other_parameters.enzyme_mass))
    network.enzyme_mass = other_parameters.enzyme_mass;
  end

end
