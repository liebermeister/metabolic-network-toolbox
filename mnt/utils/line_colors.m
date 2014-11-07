function line_colors(xxx,color_map_name)

% line_colors(xxx,color_map_name)
% 
% xxx handle to lines  (from xxx=plot(...))
% color_map_name: name of color map (default 'jet');

if ~exist('color_map_name','var'), color_map_name = 'jet'; end

if size(xxx,1)==1, xxx=xxx'; end

if length(xxx),
  eval(['farbe=' color_map_name '(' num2str(length(xxx)) ');']);
  for k=1:length(xxx)
    set(xxx(k,:),'Color',farbe(k,:))
  end
end
