function my_xticklabel(x,y,s,siz)

%my_xticklabel(x,y,s,siz)

eval(default('s','get(gca,''XTickLabel'')'));
eval(default('y','[]'));
eval(default('siz','get(gca,''Fontsize'')'));
eval(default('x','get(gca,''XTick'')'));%

if isempty(x), x = (1:length(s))'; end
if isempty(y), a = axis; y=a(3); end

set(gca,'XTick',[]);
for it=1:length(x),
   text(x(it),y,s{it},'Rotation',90,'Fontsize',siz);
end