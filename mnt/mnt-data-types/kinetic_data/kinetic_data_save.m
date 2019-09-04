function M = kinetic_data_save(kinetic_data, network, filename, options_style)

% M = kinetic_data_save(kinetic_data,network,filename,options_style)

eval(default('filename','[]','options_style','0'));

T = kinetic_data_print(kinetic_data,network,0);

M = T{1}(1,:);

for it = 1:length(T),
  if size(T{it},1)>1,
    M = [M; T{it}(2:end,:)];
  end
end

M = [[{'!!SBtab TableType="Quantity"'}, repmat({''},1,size(M,2)-1)]; M];

if length(filename),
  display(sprintf('Saving kinetic data to file %s ',filename));
  mytable(M,options_style,filename);
else,
  mytable(M,options_style)
end
