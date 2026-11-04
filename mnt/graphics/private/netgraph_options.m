function gp = netgraph_options(args, network)

% gp = netgraph_options(args, network)
%
% Graphics options given to a function in the folder 'graphics', as a struct
%
% args     cell array (e.g. varargin) containing one of
%            {}                            no options
%            {[]}                          no options
%            {gp}                          struct of options
%            {'name1', value1, ...}        name-value pairs
%            {gp, 'name1', value1, ...}    struct, followed by name-value pairs (these override the struct)
% network  (optional) network; fields of network.graphics_par count as known options
%
% Option names are checked:
%  - names that differ from a known option name only in upper and lower case are corrected
%    (e.g. 'fontsize' -> 'FontSize'); some other names are accepted as aliases (e.g. 'Rotation' -> 'fontangle')
%  - unknown names give a warning (once per name and Matlab session); the option is ignored by
%    the graphics functions. To switch the warning off: warning('off','mnt:netgraph:unknownOption')

if nargin < 2, network = struct; end

if numel(args) == 1 && iscell(args{1})
  args = args{1};
end

gp = struct;
k  = 1;
if ~isempty(args)
  if isstruct(args{1})
    gp = args{1};
    k  = 2;
  elseif isempty(args{1}) && ~ischar(args{1})
    k  = 2;
  end
end

pairs = args(k:end);
if mod(numel(pairs),2)
  error('mnt:netgraph:options', 'Graphics options must be given as a struct or as name-value pairs');
end
for it = 1:2:numel(pairs)
  name = pairs{it};
  if isstring(name), name = char(name); end
  if ~ischar(name)
    error('mnt:netgraph:options', 'Graphics options must be given as a struct or as name-value pairs');
  end
  gp.(name) = pairs{it+1};
end

if numel(gp) ~= 1
  % struct arrays (e.g. from struct('omitreactions', {...})) are left to the calling function
  return
end

gp = check_option_names(gp, network);


% ------------------------------------------------------------------

function gp = check_option_names(gp, network)

persistent warned

[names, aliases] = netgraph_option_names;

known = names;
if isstruct(network) && isfield(network,'graphics_par') && isstruct(network.graphics_par)
  known = [known, fieldnames(network.graphics_par)'];
end

fn = fieldnames(gp);
unknown = fn(~ismember(fn, known));

for it = 1:length(unknown)
  name = unknown{it};
  canonical = '';
  ind = find(strcmp(name, aliases(:,1)), 1);
  if ~isempty(ind)
    canonical = aliases{ind,2};
  else
    ind = find(strcmpi(name, names), 1);
    if ~isempty(ind), canonical = names{ind}; end
  end
  if ~isempty(canonical)
    % an option given under its correct name takes precedence
    if ~isfield(gp, canonical)
      gp.(canonical) = gp.(name);
    end
    gp = rmfield(gp, name);
  else
    if isempty(warned), warned = {}; end
    if ~ismember(name, warned)
      warned{end+1} = name;
      warning('mnt:netgraph:unknownOption', ...
              'Unknown graphics option ''%s'' (ignored). For a list of options, see ''help netgraph_draw''', name);
    end
  end
end
