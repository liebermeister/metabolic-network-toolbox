function result = load_any_table(filename,delimiter)

% list = load_any_table(filename)
%
% loads a tab-delimited file and puts it into a cell array
% filename: name of tab-delimited file containing strings and numbers

if ~exist('delimiter','var'), delimiter = sprintf('\t'); end

try

  A = {};
  fid         = fopen(filename);
  line_number = 1; 
  while 1,
    tline  = fgetl(fid);
    if ~ischar(tline), break, end
    delpos = findstr(tline,delimiter);
    delpos = [0, delpos, length(tline)+1];
    for it = 1:length(delpos)-1,
      result{line_number,it} = tline(delpos(it)+1:delpos(it+1)-1);
    end
    line_number = line_number + 1;
  end
  fclose(fid);

  result = result(:,find(sum(cellfun('length',result),1)>0));

catch 
  warning(sprintf('File %s not found',filename));
  result = [];
end