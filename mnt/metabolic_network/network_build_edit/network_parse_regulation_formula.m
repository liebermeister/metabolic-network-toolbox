function [signs, compounds] = network_parse_regulation_formula(formula)

% NETWORK_PARSE_REGULATION_FORMULA Parse a regulation formula, e.g. '+ ADP - PEP'
%
% [signs, compounds] = network_parse_regulation_formula(formula)
%
% Regulation formulae are used in the SBtab columns !KineticRegulation (Reaction
% tables: compounds that regulate the enzyme) and !ActivityCondition (Regulator
% tables: compounds that activate or inhibit the regulator). Each compound is
% preceded by a sign and a whitespace: '+' (activation), '-' (inhibition), '+-'
% (dual effect), or '?' (unknown effect; older files: no sign). The older notation '-PEP|+ADP'
% (compounds separated by '|', signs directly before the names) is also accepted.
%
% Outputs
%   signs      row vector: +1 (activation), -1 (inhibition), 0 (dual or unknown)
%   compounds  column cell array of compound names (same order as 'signs')

signs     = [];
compounds = {};
if ~ischar(formula), return; end
tokens  = strsplit(strtrim(strrep(formula,'|',' ')));
my_sign = '';
for it = 1:length(tokens)
  token = tokens{it};
  if isempty(token), continue; end
  if any(strcmp(token,{'+','-','+-','?'}))
    my_sign = token;
    continue
  end
  if strncmp(token,'+-',2) && length(token) > 2
    my_sign = '+-'; token = token(3:end);
  elseif any(token(1) == '+-?') && length(token) > 1
    my_sign = token(1); token = token(2:end);
  end
  switch my_sign
    case '+', signs(1,end+1) = 1;
    case '-', signs(1,end+1) = -1;
    otherwise, signs(1,end+1) = 0;
  end
  compounds{end+1,1} = token;
  my_sign = '';
end
