function quantity_info = data_integration_load_quantity_info(omit_quantities, quantity_info_filename)

% quantity_info = data_integration_load_quantity_info(omit_quantities,quantity_info_filename)

if ~exist('sbtab_version','file'),
  error('For this function, the SBtab toolbox must be installed');
end

eval(default('omit_quantities','[]','quantity_info_filename','[]'));

if isempty(quantity_info_filename),
  data_integration_dir   = [fileparts(which(mfilename))];
  quantity_info_filename = [ data_integration_dir '/data/quantity_info.csv'];
end

quantity_info_sbtab = sbtab_table_load(quantity_info_filename);
quantity_info       = sbtab_table_get_all_columns(quantity_info_sbtab);

keep = 1:length(quantity_info.QuantityType);
if length(omit_quantities),
  keep = setdiff(keep,label_names(omit_quantities,quantity_info.QuantityType));
end

fn = fieldnames(quantity_info);
for it = 1:length(fn),
  quantity_info.(fn{it}) =   quantity_info.(fn{it})(keep);
end

for it = 1:length(quantity_info.Symbol),
  quantity_info.symbol_index.(quantity_info.Symbol{it})=it;
end