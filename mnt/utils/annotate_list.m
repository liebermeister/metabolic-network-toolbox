function list2 = annotate_list(list1,filename)

%if key occurs several times in file,
%choose last occurrence

fid = fopen(filename);
A   = textscan(fid,'%s%s\n','delimiter','\t');
fclose(fid);

name1     = A{1};
name2     = A{2};

indices = label_names(list1,name1,'multiple');
for it =1:length(indices), 
  if     isempty(indices{it}), indices{it} = 0; 
  else, indices{it} = indices{it}(end); 
  end; 
end
indices = cell2mat(indices);
  
list2 = repmat({''},length(indices),1);

found = find(indices);
list2(found) = name2(indices(found));