function [network,N,metabolites] = network_build_from_sum_formulae(filename_reactions,filename_compounds,columns,filename_regulation)

% network = network_build_from_sum_formulae(filename_reactions,filename_compounds,columns,filename_regulation)
%
% Build matlab network structure from reactions contained in SBtab file
% (to build a network structure directly from a reaction list, use network_build_from_sum_formulae_list
%  to convert a network structure into reaction formulae, use network_print_formulae)
%
% Instead of SBtab filenames, SBtab data structures can be directly provided as function arguments
%
% attention: the syntax of the sum formulas is very strict
% example: 'A + 2 B <=> C'
% the metabolite names must not contain any spaces; double spaces or missing
% spaces will not be recognised
%
% to import the data directly from a matlab structure, the arguments
% 'filename_reactions' 'filename_compounds' need to be empty
% and the formulas have to be given in 'columns.ReactionFormula';
% other network entries can also be given
% -- OR: use the wrapper function network_build_from_sum_formulae_list --

eval(default('filename_compounds','[]','columns','[]','filename_regulation','[]'));

if length(filename_reactions)
  if ischar(filename_reactions)
    reaction_table = sbtab_table_load(filename_reactions);
  else
    reaction_table = filename_reactions;
  end
  reaction_table = sbtab_table_remove_comment_lines(reaction_table);
  columns = sbtab_table_get_all_columns(reaction_table);
end

if length(filename_reactions) * length(filename_compounds)
  if ischar(filename_compounds)
    compound_table = sbtab_table_load(filename_compounds);
  else
    compound_table = filename_compounds;
  end
  compound_table   = sbtab_table_remove_comment_lines(compound_table);
  compound_columns = sbtab_table_get_all_columns(compound_table);
end

metab_collect = {};

if exist('compound_table','var')
  compound_table = sbtab_table_remove_comment_lines(compound_table);
end
if exist('reaction_table','var')
  reaction_table = sbtab_table_remove_comment_lines(reaction_table);
end

for it = 1:length(columns.ReactionFormula)
  sum_formula              = columns.ReactionFormula{it};
  if sum(strfind(sum_formula,','))
    error(sprintf('Malformed formula "%s"',sum_formula));
  end
  if sum(strfind(sum_formula,'-'))
    warning(sprintf('Malformed formula "%s" (hyphens found)',sum_formula));
  end
  pos                      = strfind(sum_formula,'<=>');
  substrate_side           = [sum_formula(1:pos-2) ' + '];
  product_side             = [sum_formula(pos+4:end) ' + '];
  [sstoich{it},smetab{it}] = analyse_one_side(substrate_side);
  [pstoich{it},pmetab{it}] = analyse_one_side(product_side);
  metab_collect            = [metab_collect; column(smetab{it}); column(pmetab{it})];
end

% Kinetic regulation of enzymes (column KineticRegulation; older names of the
% column: MetabolicRegulation, Regulator), e.g. '+ ADP - PEP' or '-PEP|+ADP'

regulation_column = '';
for my_column = {'KineticRegulation', 'MetabolicRegulation', 'Regulator'}
  if isfield(columns, my_column{1})
    regulation_column = my_column{1};
    break
  end
end

if length(regulation_column)
  for it = 1:length(columns.(regulation_column))
    [rstoich{it,1},rmetab{it,1}] = network_parse_regulation_formula(columns.(regulation_column){it});
  end
end

if exist('compound_columns','var')
  metabolites = compound_columns.ID;
else
  metabolites = unique(metab_collect);
end

N = zeros(length(metabolites),length(columns.ReactionFormula));
for it = 1:length(columns.ReactionFormula)
  ls = label_names(smetab{it},metabolites);
  lp = label_names(pmetab{it},metabolites);
  if sum(ls==0) * length(smetab{it})
    mytable(smetab{it}')
    error(sprintf('Unknown substrate'));
  end
  if sum(lp==0) * length(pmetab{it})
    mytable(pmetab{it}')
    error(sprintf('Unknown product'));
  end
  N(ls,it) = - sstoich{it};
  N(lp,it) = pstoich{it};
end

regulation_matrix = sparse(zeros(length(columns.ReactionFormula),length(metabolites)));
kinetic_regulation_found = length(regulation_column) && any(~cellfun(@isempty, rmetab));

if length(regulation_column)
  for it = 1:length(columns.ReactionFormula)
    if isempty(rmetab{it}), continue; end
    l = label_names(rmetab{it},metabolites);
    if any(l==0)
      error('Unknown compound(s) %s in column %s (reaction %d: "%s")', strjoin(rmetab{it}(l==0)',', '), regulation_column, it, columns.(regulation_column){it});
    end
    if any(rstoich{it}==0)
      warning('Column %s, reaction %d ("%s"): regulators with dual or unspecified effect are ignored', regulation_column, it, columns.(regulation_column){it});
    end
    regulation_matrix(it,l(rstoich{it}~=0)) = rstoich{it}(rstoich{it}~=0);
  end
end

if exist('compound_columns','var')
  ll = label_names(metabolites,compound_columns.ID);
  if isfield(compound_columns,'External')
    external_ind = find(cell_string2num(compound_columns.External(ll)));
  else
    external_ind = [];
  end
else
  external_ind = [];
end

if isfield(columns,'ID')
  actions = columns.ID;
else
  if isfield(columns,'Reaction')  % for compatibility with old
    % SBtab files
    actions = columns.Reaction;
  else
    actions = numbered_names('R',length(columns.ReactionFormula));
  end
end

if isfield(columns,'IsReversible')
  reversible = sbtab_table_get_column(reaction_table,'IsReversible',1);
else
  reversible = ones(length(actions),1);
end

network = network_construct(N,reversible,external_ind,metabolites,actions,0,regulation_matrix);

if exist('compound_columns','var')

  if isfield(compound_columns,'SBML_species_ID')
    network.sbml_id_species  = compound_columns.SBML_species_ID(ll);
    compound_columns = rmfield(compound_columns,'SBML_species_ID');
  end

  if isfield(compound_columns,'SBML__species__ID')  % OLD SBTAB COMPATIBILITY
    network.sbml_id_species  = compound_columns.SBML__species__ID(ll);
    compound_columns = rmfield(compound_columns,'SBML__species__ID');
  end

  if isfield(compound_columns,'Name')
    network.metabolite_names  = compound_columns.Name(ll);
    compound_columns = rmfield(compound_columns,'Name');
  end

  if isfield(compound_columns,'NameForPlots')
    network.metabolite_NameForPlots  = compound_columns.NameForPlots(ll);
    compound_columns = rmfield(compound_columns,'NameForPlots');
  end

  if isfield(compound_columns,'Identifiers_kegg_compound')
    network.metabolite_KEGGID = compound_columns.Identifiers_kegg_compound(ll);
    compound_columns = rmfield(compound_columns,'Identifiers_kegg_compound');
  end

  if isfield(compound_columns,'Identifiers_bigg_compound')
    network.metabolite_BIGGID = compound_columns.Identifiers_bigg_compound(ll);
    compound_columns = rmfield(compound_columns,'Identifiers_bigg_compound');
  end

  if isfield(compound_columns,'Identifiers')
    dum = compound_columns.Identifiers(ll);
    for it = 1:length(dum)
      if strcmp(dum{it}(1:5),'kegg:')
        network.metabolite_KEGGID{it,1} = dum{it}(6:end);
      end
    end
    compound_columns = rmfield(compound_columns,'Identifiers');
  end

  if isfield(compound_columns,'Identifiers')
    dum = compound_columns.Identifiers(ll);
    for it = 1:length(dum)
      if strcmp(dum{it}(1:5),'bigg:')
        network.metabolite_BIGGID{it,1} = dum{it}(6:end);
      end
    end
    compound_columns = rmfield(compound_columns,'Identifiers');
  end

  if isfield(compound_columns,'IsCofactor')
    network.is_cofactor = cell_string2num(compound_columns.IsCofactor(ll));
    compound_columns = rmfield(compound_columns,'IsCofactor');
  end

end

if isfield(columns,'Name')
  network.reaction_names  = columns.Name;
  columns = rmfield(columns,'Name');
end

if isfield(columns,'NameForPlots')
  network.reaction_NameForPlots  = columns.NameForPlots;
  columns = rmfield(columns,'NameForPlots');
end

if isfield(columns,'Gene')
  network.genes  = columns.Gene;
  columns = rmfield(columns,'Gene');
end

if isfield(columns,'SBML__reaction__ID')
  network.sbml_id_reaction  = columns.SBML__reaction__ID;
  columns = rmfield(columns,'SBML__reaction__ID');
end

if isfield(columns,'Identifiers_kegg_reaction')
  network.reaction_KEGGID = columns.Identifiers_kegg_reaction;
  columns = rmfield(columns,'Identifiers_kegg_reaction');
end

if isfield(columns,'Identifiers_bigg_reaction')
  network.reaction_BIGGID = columns.Identifiers_bigg_reaction;
  columns = rmfield(columns,'Identifiers_bigg_reaction');
end

if isfield(columns,'Identifiers')
  dum = columns.Identifiers;
  for it = 1:length(dum)
    if strcmp(dum{it}(1:5),'kegg:')
      network.reaction_KEGGID{it,1} = dum{it}(6:end);
    end
  end
  columns = rmfield(columns,'Identifiers');
end

if isfield(columns,'Identifiers')
  dum = columns.Identifiers;
  for it = 1:length(dum)
    if strcmp(dum{it}(1:5),'bigg:')
      network.reaction_BIGGID{it,1} = dum{it}(6:end);
    end
  end
  columns = rmfield(columns,'Identifiers');
end

if isfield(columns,'Identifiers_ec_code')
  network.EC = columns.Identifiers_ec_code;
  columns = rmfield(columns,'Identifiers_ec_code');
end

network.formulae = columns.ReactionFormula;

if exist('compound_columns','var')
  network = join_struct(network,compound_columns);
end

network = join_struct(network,columns);


% -------------------------------------------
% Regulation (optional)

if ischar(filename_regulation) && length(filename_regulation)
  if exist(filename_regulation,'file')
    filename_regulation = sbtab_table_load(filename_regulation);
  else
    warning('Regulation file %s not found; ignoring it', filename_regulation);
    filename_regulation = [];
  end
end

if length(filename_regulation)
  % Regulation table (SBtab 1.2, TableType 'Regulation'): its regulation arrows are
  % added to the kinetic regulation from the Reaction column KineticRegulation
  [W_table, n_regulation_rows] = regulation_table_to_matrix(filename_regulation, network.metabolites, network.actions);
  W = network.regulation_matrix;
  if kinetic_regulation_found && n_regulation_rows > 0
    warning('Kinetic regulation is given both in the Reaction column %s and in a Regulation table; the information from both is combined', regulation_column);
  end
  [ind_r, ind_m] = find(W_table);
  for it = 1:length(ind_r)
    w_old = W(ind_r(it),ind_m(it));
    w_new = W_table(ind_r(it),ind_m(it));
    if w_old ~= 0 && w_old ~= w_new
      warning('Regulation of reaction %s by %s: the Regulation table and the column %s disagree; the Regulation table is used', network.actions{ind_r(it)}, network.metabolites{ind_m(it)}, regulation_column);
    end
    W(ind_r(it),ind_m(it)) = w_new;
  end
  network.regulation_matrix = W;
end

% -------------------------------------

function [sstoic,smetab] = analyse_one_side(substrate_side)

  spos = strfind(substrate_side, ' + ');
  it2 = 1;
  clear sterm
  while length(spos)
    try
      sterm{it2} = substrate_side(1:spos(1)-1);
      wpos = strfind(sterm{it2},' ');
      if wpos
        sstoic(it2) = eval(sterm{it2}(1:wpos(1)-1));
        smetab{it2} = sterm{it2}(wpos(1)+1:end);
      else
        sstoic(it2) = 1;
        smetab{it2} = sterm{it2};
      end
      substrate_side = substrate_side(spos(1)+3:end);
      spos = strfind(substrate_side, ' + ');
      it2 = it2+1;
    catch
      substrate_side(1:spos(1)-1)
      error('Parsing error');
    end
  end

  if length(smetab{1})==0
    smetab=[];
  end


% -------------------------------------

function [W, n] = regulation_table_to_matrix(regulation_table, metabolites, reactions)

% Regulation matrix (reactions x compounds; +1 activation, -1 inhibition) from an
% SBtab Regulation table with columns
%   !Regulator      (older name: !Compound)    regulating compound
%   !Target         (older names: !TargetReaction, !Reaction)   regulated reaction
%   !RegulatorType, !TargetType  (optional) entity types; rows with types other than
%                   'metabolite' (regulator) and 'enzyme' (target), e.g. regulation of
%                   gene expression, are skipped
%   !RegulationType          'activation', 'inhibition', or 'dual'
%   !RegulationType:SBO      SBO term, e.g. SBO:0000459 (stimulator), SBO:0000020 (inhibitor)
%   !RegulatoryCoefficient   signed Hill coefficient (+ activation, - inhibition; NaN: unknown or dual)
% The sign is taken from the first of these three columns that has an entry;
% if the columns disagree, a warning is given. Regulators with dual or unspecified
% effect, and unknown reactions or compounds, are ignored with a warning.
% (The magnitudes of regulatory coefficients are not stored yet.)

nr = length(reactions);
nm = length(metabolites);
W  = sparse(nr,nm);
n  = 0;

reaction_column  = first_existing_column(regulation_table, {'Target','TargetReaction','Reaction'});
regulator_column = first_existing_column(regulation_table, {'Regulator','Compound'});
if isempty(reaction_column) || isempty(regulator_column)
  warning('Regulation table without columns Target and Regulator; ignoring it');
  return
end
my_reactions  = sbtab_table_get_column(regulation_table, reaction_column);
my_regulators = sbtab_table_get_column(regulation_table, regulator_column);
n = length(my_reactions);

sign_type = nan(n,1);
sign_sbo  = nan(n,1);
sign_coef = nan(n,1);
if sbtab_table_has_column(regulation_table,'RegulationType')
  my_types = sbtab_table_get_column(regulation_table,'RegulationType');
  for it = 1:n
    sign_type(it) = regulation_type_sign(my_types{it});
  end
end
if sbtab_table_has_column(regulation_table,'RegulationType:SBO')
  my_sbo = sbtab_table_get_column(regulation_table,'RegulationType:SBO');
  for it = 1:n
    sign_sbo(it) = regulation_type_sign(my_sbo{it});
  end
end
if sbtab_table_has_column(regulation_table,'RegulatoryCoefficient')
  my_coef = sbtab_table_get_column(regulation_table,'RegulatoryCoefficient',1);
  sign_coef = sign(my_coef(:));
  sign_coef(sign_coef == 0) = nan;
end

ind_reactions  = label_names(my_reactions, reactions);
ind_regulators = label_names(my_regulators, metabolites);

% rows that describe no kinetic regulation of an enzyme by a compound
is_kinetic = true(n,1);
type_columns = {'RegulatorType', 'TargetType'};
kinetic_types = {'metabolite', 'enzyme'};
for ic = 1:2
  if sbtab_table_has_column(regulation_table, type_columns{ic})
    my_entity_types = sbtab_table_get_column(regulation_table, type_columns{ic});
    for it = 1:n
      t = my_entity_types{it};
      if ischar(t) && ~isempty(strtrim(t)) && ~strcmpi(strtrim(t), kinetic_types{ic})
        is_kinetic(it) = false;
      end
    end
  end
end

for it = 1:n
  if ~is_kinetic(it), continue; end
  signs = [sign_type(it), sign_sbo(it), sign_coef(it)];
  given = signs(~isnan(signs));
  if isempty(given) || given(1) == 0
    warning('Regulation table, row %d (%s, %s): dual or unspecified effect; ignored', it, my_reactions{it}, my_regulators{it});
    continue
  end
  if any(given ~= given(1))
    warning('Regulation table, row %d (%s, %s): RegulationType, RegulationType:SBO and RegulatoryCoefficient disagree; using the first given', it, my_reactions{it}, my_regulators{it});
  end
  if ind_reactions(it) == 0 || ind_regulators(it) == 0
    warning('Regulation table, row %d: unknown reaction %s or compound %s; ignored', it, my_reactions{it}, my_regulators{it});
    continue
  end
  W(ind_reactions(it), ind_regulators(it)) = given(1);
end

% -------------------------------------

function s = regulation_type_sign(entry)

% +1 (activation), -1 (inhibition), 0 (dual) or nan (empty or unknown) for an entry of
% the column RegulationType ('activation', 'inhibition', 'dual'; older files may contain
% 'activator', 'inhibitor' or entries like 'SBO:0000020 - inhibitor') or RegulationType:SBO

activator_sbo = [459, 21, 461, 462, 170, 171, 656];
inhibitor_sbo = [20, 206, 207, 536, 537, 169];
dual_sbo      = [595];

s = nan;
if ~ischar(entry) || isempty(strtrim(entry)), return; end
entry = lower(strtrim(entry));
sbo = regexp(entry, 'sbo:0*(\d+)', 'tokens', 'once');
if ~isempty(sbo)
  sbo = str2double(sbo{1});
  if any(sbo == activator_sbo), s = 1;
  elseif any(sbo == inhibitor_sbo), s = -1;
  elseif any(sbo == dual_sbo), s = 0;
  end
  return
end
switch entry
  case {'activation','activator','+'},  s = 1;
  case {'inhibition','inhibitor','-'},  s = -1;
  case {'dual','+-'},                   s = 0;
end

% -------------------------------------

function column_name = first_existing_column(table, names)

column_name = '';
for it = 1:length(names)
  if sbtab_table_has_column(table, names{it})
    column_name = names{it};
    return
  end
end
