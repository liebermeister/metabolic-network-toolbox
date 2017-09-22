function commandstring = default(varargin)

% commandstring = default(varname1,defaultvalue1,varname2,defaultvalue2,...)
%
% usage eval(default('a','1','name','''Gonzo'''));

commandstring = '';

for it = 1:length(varargin)/2,

varname = varargin{it*2-1};
value_string = varargin{it*2};
commandstring = [commandstring, 'if ~exist(''',varname,''',''var''), ',varname,' = ',value_string,'; end; '];

end
