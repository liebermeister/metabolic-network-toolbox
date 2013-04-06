function quantity_info = data_integration_load_quantity_info(omit_quantities,quantity_info_filename)

% quantity_info = data_integration_load_quantity_info(omit_quantities,quantity_info_filename)

eval(default('omit_quantities','[]','quantity_info_filename','[]'));

if isempty(quantity_info_filename),
data_integration_dir = [fileparts(which(mfilename))];
quantity_info_filename = [ data_integration_dir '/quantity_info.tsv'];
end

quantity_info_sbtab = sbtab_load_table(quantity_info_filename);
quantity_info       = quantity_info_sbtab.column.column;

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