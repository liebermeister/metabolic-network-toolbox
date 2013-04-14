function x_centre = find_polytope_centre(M_eq,b_eq,M_ineq,b_ineq,x_lower,x_upper, x_try)

% x_centre = find_polytope_centre(M_eq, b_eq, M_ineq, b_ineq, x_lower, x_upper, x_try)
%
% x_try >= x_lower
% x_try <= x_upper
% M_ineq * x_try <= b_ineq
% M_eq * x_try == b_eq

if 0,
 %% test
 M_eq    = [];
 b_eq    = [];
 M_ineq  = [];
 b_ineq  = [];
 x_lower = -1;
 x_upper =  1;
 x = find_polytope_centre(M_eq,b_eq,M_ineq,b_ineq,x_lower,x_upper);
end

x_length = length(x_lower);

if x_length == 0, 
  x_centre = []; 
  warning('x vector has zero lenght');
  return
end

% using this lead to problems before (no feasible solution found)!!
% opt = optimset('MaxIter',10^10,'LargeScale','off','Display','off');

opt = optimset('Display','off');

for it = 1:x_length,
  z = zeros(x_length,1);
  z(it) = 1;
  [solution1,value,exitflag1] = linprog(z,M_ineq,b_ineq,M_eq,b_eq,x_lower,x_upper,[],opt);
  z(it) = -1;
  [solution2,value,exitflag2] = linprog(z,M_ineq,b_ineq,M_eq,b_eq,x_lower,x_upper,[],opt);
  if [exitflag1 ~=1] + [exitflag2 ~=1],
    [exitflag1  exitflag2]
    error('No feasible point found');
  else
   x1(:,it) = solution1; 
   x2(:,it) = solution2; 
  end
end

x_centre = mean([x1 x2],2);
