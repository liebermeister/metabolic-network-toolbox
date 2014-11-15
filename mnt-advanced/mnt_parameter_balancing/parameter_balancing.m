function [res, exitflag] = parameter_balancing(task, quantity_info, options)

% [res, exitflag] = parameter_balancing(task, quantity_info, options)
%
% Parameter balancing based on precalculated vectors and matrices 
% in structure 'task' (built by function parameter_balancing_task.m)
%
% quantity_info: (optional) table of biochemical quantities 
% (from data_integration_load_quantity_info.m)

eval(default('quantity_info','[]','options','struct'));

if isempty(quantity_info),
  quantity_info = data_integration_load_quantity_info;
end

options_default.insert_pseudo_values = 0;

options = join_struct(options_default,options);

q_prior.cov_inv = diag(sparse(1./task.q.prior.std.^2));

if find(task.xdata.std==0),
  warning('Zero standard deviations found. Replcaing them by 10^-10.'); 
  task.xdata.std(task.xdata.std==0) = 10^-5;
end

xdata.cov_inv   = diag(sparse(1./task.xdata.std.^2));

xall_pseudo.cov_inv = diag(sparse(1./task.xall.pseudo.std.^2));

%flag_bounds = length(task.xlower.value) + length(task.xupper.value) > 0;
flag_bounds = 1;

if length(task.Q_xdata_q),
  q_posterior_cov_inv = q_prior.cov_inv + [ task.Q_xdata_q' * xdata.cov_inv * task.Q_xdata_q ];
  TT                  = task.Q_xdata_q' * xdata.cov_inv * task.xdata.mean + q_prior.cov_inv * task.q.prior.mean;
else,
  q_posterior_cov_inv = q_prior.cov_inv;
  TT                  = q_prior.cov_inv * task.q.prior.mean;
end

if options.insert_pseudo_values, 
  q_posterior_cov_inv = q_posterior_cov_inv + task.Q_xall_q' * xall_pseudo.cov_inv * task.Q_xall_q;
  TT                  = TT + task.Q_xall_q' * xall_pseudo.cov_inv * task.xall.pseudo.mean;
end
q_posterior.mean    = q_posterior_cov_inv \ TT;

if flag_bounds, 
  xconstraints = [-task.xlower.value_nat; task.xupper.value_nat;];
  Qconstraints = [-task.Q_xlower_q;   task.Q_xupper_q;];
  if options.insert_pseudo_values, 
    xconstraints = [xconstraints; -task.xall.lower_nat; task.xall.upper_nat;];
    Qconstraints = [Qconstraints; -task.Q_xall_q;       task.Q_xall_q];
  end
  %% [Qconstraints * q_posterior.mean, xconstraints]
  epsilon = 10^-10;
  opt = optimset('Display','off','Algorithm','active-set','MaxIter',10^8);
  [q_posterior.mode,fval,exitflag] = quadprog(full(q_posterior_cov_inv), full(-q_posterior_cov_inv * q_posterior.mean), full(Qconstraints), xconstraints - epsilon,[],[],[],[],[],opt);
  %% [Qconstraints * q_posterior.mode, xconstraints]

else
  q_posterior.mode = q_posterior.mean;
end

if exitflag~=1, 
  exitflag
  if exitflag == 0, 
    error('Problem in quadratic programming; increase iteration number');
  end
  error('Problem in quadratic programming; possible reasons are too tight constraints (for concentrations, reaction affinities etc, or an infeasible flux distribution');
  
end

q_posterior.cov       = inv(q_posterior_cov_inv);
q_posterior.std       = sqrt(diag(q_posterior.cov));

xmodel_posterior.mode = task.Q_xmodel_q * q_posterior.mode;
xmodel_posterior.mean = task.Q_xmodel_q * q_posterior.mean;
xmodel_posterior.cov  = task.Q_xmodel_q * q_posterior.cov * task.Q_xmodel_q';
xmodel_posterior.std  = sqrt(diag(xmodel_posterior.cov));

if length(task.xdata.mean),
  xdata_posterior.mode  = task.Q_xdata_q * q_posterior.mode;
  xdata_posterior.mean  = task.Q_xdata_q * q_posterior.mean;
  xdata_posterior.cov   = task.Q_xdata_q * q_posterior.cov * task.Q_xdata_q';
  xdata_posterior.std   = sqrt(diag(xdata_posterior.cov));
end

res.q_posterior.mode      = q_posterior.mode     ;
res.q_posterior.mean      = q_posterior.mean     ;
res.q_posterior.cov       = q_posterior.cov      ;
res.q_posterior.std       = q_posterior.std      ;
res.xmodel_posterior.mode = xmodel_posterior.mode;
res.xmodel_posterior.mean = xmodel_posterior.mean;
res.xmodel_posterior.cov  = xmodel_posterior.cov ;
res.xmodel_posterior.std  = xmodel_posterior.std ;

if length(task.xdata.mean),
  res.xdata_posterior.mode  = xdata_posterior.mode ;
  res.xdata_posterior.mean  = xdata_posterior.mean ;
  res.xdata_posterior.cov   = xdata_posterior.cov  ;
  res.xdata_posterior.std   = xdata_posterior.std  ;
end

% -----------------------------------------------------------------------
% translate the posterior mean values back to kinetic constants

[nr,nm,nx,ind_KM,ind_KA,ind_KI,nKM,nKA,nKI] = network_numbers(task.network);

num_model = parameter_balancing_quantity_numbers(task.model_quantities,quantity_info,task.network);

for it = 1:length(task.model_quantities),
  my_quantity     = task.model_quantities{it};
  ind             = find(strcmp(my_quantity,quantity_info.QuantityType));
  my_scaling      = quantity_info.Scaling{ind};
  my_symbol       = quantity_info.Symbol{ind};
  my_indices      = task.xmodel.indices.(my_symbol);
  my_x_mode       = xmodel_posterior.mode(my_indices);
  my_x_mean       = xmodel_posterior.mean(my_indices);
  my_x_std        = xmodel_posterior.std(my_indices);

  if strcmp( my_scaling, 'Logarithmic'),
     my_x_mode            = exp(my_x_mode);
     my_x_median          = exp(my_x_mean);
    [my_x_mean, my_x_std] = lognormal_log2normal(my_x_mean, my_x_std);
  else
    my_x_median          = my_x_mean;
  end
  
  switch my_symbol, 
    case 'KM', indices = ind_KM;
    case 'KA', indices = ind_KA;
    case 'KI', indices = ind_KI;
  end
  switch my_symbol, 
    case {'KM', 'KA', 'KI'},
      dum = sparse(zeros(nr,nm));
      dum_mode   = dum; dum_mode(indices)   = my_x_mode;   my_x_mode   = dum_mode;
      dum_mean   = dum; dum_mean(indices)   = my_x_mean;   my_x_mean   = dum_mean;
      dum_std    = dum; dum_std(indices)    = my_x_std;    my_x_std    = dum_std;
      dum_median = dum; dum_median(indices) = my_x_median; my_x_median = dum_median;
  end

  kpm_mode.(my_symbol)   = my_x_mode;
  kpm_mean.(my_symbol)   = my_x_mean;
  kpm_std.(my_symbol)    = my_x_std;
  kpm_median.(my_symbol) = my_x_median;
end

% -------------------------------------------------------------

res.kinetics_posterior_mode   = kpm_mode;
res.kinetics_posterior_median = kpm_median;
res.kinetics_posterior_mean   = kpm_mean;
res.kinetics_posterior_std    = kpm_std;

