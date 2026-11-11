function [keq, info] = mnt_keq_from_kinetic_data(network, kinetic_data_file, options)

% [keq, info] = mnt_keq_from_kinetic_data(network, kinetic_data_file, options)
%
% Equilibrium constants for a model, obtained by parameter balancing from kinetic and
% thermodynamic data (requires the parameter balancing toolbox)
%
% Inputs
%   network            network structure
%   kinetic_data_file  SBtab file with kinetic/thermodynamic data (table 'ParameterData' or similar,
%                      as read by parameter_balancing)
%   options            struct with (optional) fields
%     .model_file         model passed to parameter_balancing: SBtab or SBML file (default: the network struct)
%     .model_name         (default 'model')
%     .prior_file         prior table for parameter balancing (default: mnt_prior_file)
%     .output_file        file for the balanced parameters (default []: not saved; parameter balancing
%                         also writes files with the endings _preprocessed_data.tsv, _options.tsv)
%     .match_data_by      'KeggId' (default), 'ModelElementId' or 'BiggId' (data rows -> model reactions)
%     .orient_data        (default 1) convert the data to the reaction orientation of the model first
%                         (sbtab_orient_data_files); see mnt_map_state_data
%     .model_orientation  (default []) name of the reaction orientation of the model
%     .verbose            (default 0)
%
% Outputs
%   keq   vector of equilibrium constants (median values of the balanced parameters), one per reaction
%   info  struct with fields .balanced_parameters (SBtab table), .kinetic_data (balanced parameters, as
%         returned by kinetic_data_load), .orientation_report
%
% The standard concentration of the equilibrium constants is the one of the data and prior files.

eval(default('options','struct'));

options_default = struct('model_file', [], 'model_name', 'model', 'prior_file', [], 'output_file', [], ...
                         'match_data_by', 'KeggId', 'orient_data', 1, 'model_orientation', [], 'verbose', 0);
options = join_struct(options_default, options);

if isempty(options.prior_file), options.prior_file = mnt_prior_file; end
model = options.model_file;
if isempty(model), model = network; end

info.orientation_report = [];
if options.orient_data
  [kinetic_data_file, info.orientation_report] = sbtab_orient_data_files(network, kinetic_data_file, ...
       struct('match_data_by', options.match_data_by, 'model_orientation', options.model_orientation, 'verbose', options.verbose));
end

pb_options                          = parameter_balancing_options;
pb_options.run_variability_analysis = 0;
pb_options.use_sbml_ids             = 0;
pb_options.use_kegg_ids             = strcmp(options.match_data_by, 'KeggId');
pb_options.flag_check               = options.verbose;
pb_options.verbose                  = options.verbose;

balanced_parameters = parameter_balancing(model, options.output_file, kinetic_data_file, options.prior_file, [], options.model_name, pb_options);

if ~isempty(options.output_file)
  sbtab_table_save(balanced_parameters, options.output_file);
end

parameter_prior   = parameter_balancing_prior([], options.prior_file, 0);
info.kinetic_data = kinetic_data_load([], parameter_prior, network, balanced_parameters, struct('enforce_ranges', 0));
info.balanced_parameters = balanced_parameters;

keq = info.kinetic_data.Keq.median;
if any(isnan(keq))
  warning('Equilibrium constants missing for %d reactions', sum(isnan(keq)));
end
