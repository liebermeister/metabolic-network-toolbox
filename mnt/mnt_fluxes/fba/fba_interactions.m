function [f_ref,f_single,f_double,v_ref,v_single,v_double,df1,df2] = fba_interactions(network,fba_constraints,inhibition_factor)

% [f_ref,f_single,f_double,v_ref,v_single,v_double,df1,df2] = fba_interactions(network,fba_constraints,inhibition_factor)
%
% f_ref      unperturbed objective value
% f_single   single-perturbation objectives (vector)
% f_double   double-perturbation objectives (matrix)
% v_ref      unperturbed flux mode
% v_single   single-perturbation flux modes (matrix)
% v_double   double-perturbation flux modes (tensor)
% df1        single-perturbation objective deficits (vector)
% df2        double-perturbation synergistic objective deficits (matrix)

epsilon = 10^-8;

[nm,nr]        = size(network.N);
[v_ref, f_ref] = fba(network,fba_constraints);

for it1 = 1:nr,

  %% first order
  fba_constraints_1 = fba_constraints;
  if v_ref(it1)>0,
    fba_constraints_1.v_max(it1) = inhibition_factor * v_ref(it1);
  end
  if v_ref(it1)<0,
    fba_constraints_1.v_min(it1) = inhibition_factor * v_ref(it1);
  end

  %% second order
  [v_single(:,it1),f_single(it1,1)] = fba(network,fba_constraints_1);
  for it2 = 1:it1-1,
    fba_constraints_2 = fba_constraints_1;
    if v_ref(it2)>0,
      fba_constraints_2.v_max(it2) = inhibition_factor * v_ref(it2);
    end
    if v_ref(it2)<0,
      fba_constraints_2.v_min(it2) = inhibition_factor * v_ref(it2);
    end
    [v_double(:,it1,it2),f_double(it1,it2)] = fba(network,fba_constraints_2);
  end

  % diagonal elements (square perturbation!)
  it2 = it1;
  fba_constraints_2 = fba_constraints_1;
  if v_ref(it2)>0,
    fba_constraints_2.v_max(it2) = inhibition_factor^2 * v_ref(it2);
  end
  if v_ref(it2)<0,
    fba_constraints_2.v_min(it2) = inhibition_factor^2 * v_ref(it2);
  end
  [v_double(:,it1,it2),f_double(it1,it2)] = fba(network,fba_constraints_2);
end

f_double = f_double + f_double' - diag(diag(f_double));

f_single(abs(f_single)<epsilon) = 0;
f_double(abs(f_double)<epsilon) = 0;

df1 = f_single - f_ref; 
df2 = f_double - f_ref - repmat(df1,1,nr) - repmat(df1',nr,1);
