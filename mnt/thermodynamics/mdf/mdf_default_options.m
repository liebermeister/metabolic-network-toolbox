function mdf_options = mdf_default_options(network,opt)

% mdf_options = mdf_default_options(network,opt)
%
% Default options for the max-min driving force method (see mdf.m).
% The default values are read from the SBtab Config file 'max-min-driving-force_Config.tsv'
% (in the same folder as this function; TableID 'ConfigureMDF'). To change the defaults, edit this file.
%
% Options (concentrations are in mM):
%
% mdf_options.c_min                 lower bounds for metabolite concentrations (vector, nm x 1; nan if not specified)
% mdf_options.c_max                 upper bounds for metabolite concentrations (vector, nm x 1; nan if not specified)
% mdf_options.c_fix                 fixed metabolite concentrations (vector, nm x 1; nan if not specified)
% mdf_options.weight_by_fluxes      flag: during optimisation, weight the forces by flux values
% mdf_options.ind_ignore_reactions  indices of reactions to be ignored in mdf
%
% where nm = # metabolites (rows of full stoichiometric matrix N)
%
% The argument 'opt' (optional struct) can contain the fields
%   opt.c_min_default         default lower bound, used for all metabolites
%   opt.c_max_default         default upper bound, used for all metabolites
%   opt.ind_ignore_reactions  indices of reactions to be ignored
%   opt.fixed_met             fixed metabolite concentrations, given as a struct
%                             with elements opt.fixed_met.[METABOLITE_NAME] == [METABOLITE_VALUE]

eval(default('opt','struct'));

config = load_config_table(fullfile(fileparts(mfilename('fullpath')), 'max-min-driving-force_Config.tsv'));

opt_def.c_min_default = config.c_min_default; % mM
opt_def.c_max_default = config.c_max_default; % mM

opt = join_struct(opt_def,opt);

[nm,nr] = size(network.N);

mdf_options.c_min = opt.c_min_default * ones(nm,1);
mdf_options.c_max = opt.c_max_default * ones(nm,1);
mdf_options.c_fix = config.c_fix * ones(nm,1);

mdf_options.weight_by_fluxes = config.weight_by_fluxes;

mdf_options.ind_ignore_reactions = config.ind_ignore_reactions;
if isfield(opt,'ind_ignore_reactions')
  mdf_options.ind_ignore_reactions = opt.ind_ignore_reactions;
end

% set predefined concentrations

if isfield(opt,'fixed_met')
  fn = fieldnames(opt.fixed_met);
  for it = 1:length(fn)
    mdf_options.c_fix(label_names(fn{it},network.metabolites)) = opt.fixed_met.(fn{it});
  end
end
