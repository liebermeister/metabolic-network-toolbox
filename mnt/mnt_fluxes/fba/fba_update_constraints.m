function fba_constraints = fba_update_constraints(fba_constraints)

% fba_constraints = fba_update_constraints(fba_constraints)
%
% make fba constraints consistent

fba_constraints.v_min(find([fba_constraints.v_sign == 1] .*        [fba_constraints.v_min<0])) = 0;
fba_constraints.v_min(find([fba_constraints.v_sign == 1] .* ~isfinite(fba_constraints.v_min))) = 0;
fba_constraints.v_max(find([fba_constraints.v_sign ==-1] .*        [fba_constraints.v_max>0])) = 0;
fba_constraints.v_max(find([fba_constraints.v_sign ==-1] .* ~isfinite(fba_constraints.v_max))) = 0;

fba_constraints.v_sign(find(fba_constraints.v_min>0)) = 1;
fba_constraints.v_sign(find(fba_constraints.v_max<0)) = -1;

fba_constraints.v_sign(fba_constraints.v_fix==0) =  0;
fba_constraints.v_sign(fba_constraints.v_fix>0)  =  1;
fba_constraints.v_sign(fba_constraints.v_fix<0)  = -1;
fba_constraints.v_fix(fba_constraints.v_sign==0) =  0;

if sum(fba_constraints.v_fix(fba_constraints.v_sign>0) == 0) ...
      + sum(fba_constraints.v_fix(fba_constraints.v_sign<0) == 0),
  error('contradicting constraints');
end

fba_constraints.zv(isnan(fba_constraints.zv)) = 0;
