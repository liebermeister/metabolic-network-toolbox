function data = biotable_read(basename,n_name_columns,p)

% data = biotable_read(basename,n_name_columns,p)
%
% read standard table format for dynamic data and create matlab structure
%
% p.extension: file extension to be added if basename has none of the extensions
%              .txt, .csv, .tsv (default '.tsv')
%
% For preprocessing the data, see biotable_smooth_mean and biotable_omit_outliers

options_default = struct('extension','.tsv');

eval(default('p','struct'));
p = join_struct(options_default,p);

if ~endsWith(basename, {'.txt','.csv','.tsv'})
  basename = [basename, p.extension];
end

A    = load_any_table(basename);
A    = A(~strncmp(A(:,1),'#',1),:);

titles = A(1,:);
names  = A(3:end,1:n_name_columns);

data   = struct;

for it = 1:n_name_columns
  tt   = titles{it};
  if strcmp(tt(1),'!'), tt = tt(2:end); end
  data = setfield(data,tt,names(:,it));
end

data.SampleName =  titles(n_name_columns+1:2:end)';
data.SampleTime =  cell_string2num(A(2,n_name_columns+1:2:end-1))';
data.DataMean   =  cell_string2num(A(3:end,n_name_columns+1:2:end-1));
data.DataStd    =  cell_string2num(A(3:end,n_name_columns+2:2:end));

data.Info = {};
