function c = join_struct(a,b)

% c = join_struct(a,b)
%
% start with structure a and insert all fields from structure b

c = a;

ff = fieldnames(b);

for it = 1:length(ff),
  c = setfield(c,ff{it},getfield(b,ff{it}));
end