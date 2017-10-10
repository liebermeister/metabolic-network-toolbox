function parameter_prior = biochemical_parameter_prior(omit_quantities, parameter_prior_filename)

% parameter_prior = biochemical_parameter_prior(omit_quantities,parameter_prior_filename)

if ~exist('sbtab_version','file'),
  error('For this function, the SBtab toolbox must be installed');
end

eval(default('omit_quantities','[]','parameter_prior_filename','[]'));

if isempty(parameter_prior_filename),
  data_integration_dir     = [fileparts(which(mfilename))];
  parameter_prior_filename = [ data_integration_dir '/ecm_parameter_prior_10-2017.tsv'];
end

parameter_prior_sbtab = sbtab_table_load(parameter_prior_filename);
parameter_prior       = sbtab_table_get_all_columns(parameter_prior_sbtab);

keep = 1:length(parameter_prior.QuantityType);

if length(omit_quantities),
  keep = setdiff(keep,label_names(omit_quantities,parameter_prior.QuantityType));
end

fn = fieldnames(parameter_prior);
for it = 1:length(fn),
  parameter_prior.(fn{it}) = parameter_prior.(fn{it})(keep);
end

for it = 1:length(parameter_prior.Symbol),
  parameter_prior.symbol_index.(parameter_prior.Symbol{it})=it;
end