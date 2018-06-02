function N = cell_string2num(C)

% N = cell_string2num(C)

N = zeros(size(C));
for i=1:size(C,1),
  for k=1:size(C,2),
 if length( str2num(C{i,k})),   
     N(i,k)=     str2num(C{i,k});
     else, 
     N(i,k)=nan;
end     
  end
end