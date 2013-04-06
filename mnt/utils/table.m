%result = table(names,rownumbers,filename) 
%
%display cell array of strings to file
%flag rownumbers: 0 (none), 1 (yes) ,'tex' (tex style)
%(if filename is given -> output to file)

function res = table(names,rownumbers,filename)

eval(default('rownumbers','1','filename','[]'));

[nlines,nfields] = size(names);

res = '';

if nlines*nfields == 0, warning('Empty table'); end

switch rownumbers,
    
  case 0,
    for k=1:nlines
      for l=1:nfields-1
        res = [res, sprintf('%s\t',str(names{k,l}))];
      end
      res = [res, sprintf('%s\n', str(names{k,nfields}))];
    end
    
  case 'tex',
    for k=1:nlines
      for l=1:nfields-1
        res = [res, sprintf('%s \& ',str(names{k,l}))];
      end
      res = [res, sprintf('%s\\\\\n',str(names{k,nfields}))];
    end
    
  otherwise,
    for k=1:nlines
      res = [res, sprintf('%d\t',k)];
      for l=1:nfields-1
        res = [res, sprintf('%s\t',str(names{k,l}))];
      end
      res = [res, sprintf('%s\n', str(names{k,nfields}))];
    end
end
  
if ~isempty(filename),
  file = fopen(filename,'w');
  fprintf(file,'%s',res);
  fclose(file);
end

% ----------------------------------------------------------------------

function res = str(text)

if iscell(text), 
  res = text{1}; 
  for i=2:length(text), res=[res char(9) text{i}]; end
elseif isnumeric(text),
  res = num2str(text,3);
else,
  res = char(text);
end 
