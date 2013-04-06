function [v,C] = eba_make_feasible(v,N,eba_condition,C,ind_ignore,cycle_method,ind_ext)

% v = eba_make_feasible(v,N,eba_condition,C,ind_ignore)
%
% Replace a given flux distribution v by a similar, 
% but thermodynamically feasible flux distribution
%
% eba_condition: 'loose', 'strict'
% ind_ignore: indices of reactions that can be ignored in calculation of efm
%
% C is a matrix of thermodynamically unfeasible elementary cycles
%
% if it is not provided, it is computed by one of two methods
%    cycle_method: 'beard', 'efmtool'
% ind_ext: optional; if given -> projection to stationary flux distribution

eval(default('eba_condition','''strict''','C','nan','ind_ignore','[]','cycle_method','''beard'''));

v_orig = v;

if isnan(C),
  display('Computing elementary cycles. This may take a while');
  switch cycle_method,
    case 'beard',
      C = eba_my_cycles(N,C,ind_ignore);
    case 'efmtool',
      [nm,nr] = size(N);
      network = network_construct(full(N),ones(nr,1),[]);
      C = network_efmtool(network,'total',ind_ignore);
  end
end

feasible = 0;

[feasible,CC,ind_violated] = eba_feasible(sign(v),N,C,ind_ignore,eba_condition);

% adjust flux distribution by substracting a multiple of each 
% cycle involved in a violation
% repeat this until a feasible solution is reached

while ~feasible,
  fprintf('.')
  %figure(1); netgraph_concentrations(network_CoHid,[],v);
  %figure(2); netgraph_concentrations(network_CoHid,[],C(:,ind_violated(1)));
  for it = 1:length(ind_violated),
    a   = C(:,ind_violated(it)).*v;
    if sum(a~=0),
      switch eba_condition,
        case 'loose',  dum = min(abs(a(a~=0))) * sign(sum(a));
          %% choose prefactor such that the smallest flux is shifted to zero
        case 'strict', dum = mean(abs(a(a~=0))) * sign(sum(a));
          %% choose prefactor such that some fluxes change their signs
      end
      v   = v - C(:,ind_violated(it))*dum;
    end
  end
  v(abs(v)<10^-10) = 0;
  if exist('ind_ext','var'),
    v = project_fluxes(N, ind_ext, v);
  end
  %% check for feasibility
  [feasible,CC,ind_violated] = eba_feasible(sign(v),N,C,ind_ignore,eba_condition);
  ind_violated
end

if feasible, display('Feasible flux distribution found'); end
