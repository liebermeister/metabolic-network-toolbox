function x = num2cellstr(y)

x= cellstr(num2str(y));
 for it=1:length(x),
   x{it} = strrep(x{it},' ','');
end
