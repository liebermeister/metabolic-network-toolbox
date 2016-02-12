function result = parameter_balancing_process_matrix(quantity_info)

% process the formulae
ind_indep_quantities = find(strcmp(quantity_info.Dependence,'Basic'));
n_quantities         = size(quantity_info.QuantityType,1);
n_indep_quantities   = length(ind_indep_quantities);

matrix_terms = repmat({},n_quantities,n_indep_quantities);
for it = 1:n_quantities,
  terms = strsplit(',',strrep(quantity_info.MatrixInfo{it}(2:end-1),' ',''));  % omit outer brackets
   % omit inner brackets
  for itt = 1:length(terms), 
    terms(itt) = strrep(terms(itt),'[',''); 
    terms(itt) = strrep(terms(itt),']',''); 
  end
  matrix_terms(it,:) = terms;
end

% check if basic parameters are in fact independent ???

result.all_quantities         = quantity_info.QuantityType;
result.all_basic_quantities   = quantity_info.QuantityType(ind_indep_quantities);
result.matrix_terms           = matrix_terms;
