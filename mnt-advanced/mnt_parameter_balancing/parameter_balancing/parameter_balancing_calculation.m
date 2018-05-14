function [result, exitflag] = parameter_balancing_calculation(task, parameter_prior, options)
  
% [result, exitflag] = parameter_balancing(task, parameter_prior, options)
%
% Parameter balancing based on vectors and matrices in data structure 'task' 
%
% Function inputs (mandatory)
%  task:  parameter balancing task (see parameter_balancing_task.m)
%
% Function inputs (optional)
%  parameter_prior: table of biochemical quantities (see biochemical_parameter_prior.m)
%  options:       options table with fields 
%                   options.use_pseudo_values (default 0)
%                   options.n_samples            (default 0)
%                     number of random samples from posterior 
%                     for output result.kinetics_posterior_samples
%                   options.use_bounds_from_quantity_table
%                   options.fix_Keq_in_sampling


eval(default('parameter_prior','[]','options','struct'));

if isempty(parameter_prior),
  parameter_prior = biochemical_parameter_prior;
end

options = join_struct(parameter_balancing_default_options,options);

q_prior.cov_inv = diag(sparse(1./task.q.prior.std.^2));

if find(task.xdata.std==0),
  warning('Zero standard deviations found. Replacing them by 10^-10.'); 
  task.xdata.std(task.xdata.std==0) = 10^-5;
end

xdata.cov_inv   = diag(sparse(1./task.xdata.std.^2));

xall_pseudo.cov_inv = diag(sparse(1./task.xpseudo.std.^2));

flag_bounds = length(task.xlower.value_nat) + length(task.xupper.value_nat) > 0;


% ------------------------------------------------------------
% Posterior mean and covariance matrix (without constraints)

% Terms from prior
q_posterior_cov_inv = q_prior.cov_inv;
TT                  = q_prior.cov_inv * task.q.prior.mean;

%% Add pseudo values terms
if options.use_pseudo_values,
  q_posterior_cov_inv = q_posterior_cov_inv + task.Q_xpseudo_q' * xall_pseudo.cov_inv * task.Q_xpseudo_q;
  TT                  = TT                  + task.Q_xpseudo_q' * xall_pseudo.cov_inv * task.xpseudo.mean;
end

% Add data terms
if length(task.Q_xdata_q),
  q_posterior_cov_inv = q_posterior_cov_inv + task.Q_xdata_q' * xdata.cov_inv * task.Q_xdata_q;
  TT                  = TT                  + task.Q_xdata_q' * xdata.cov_inv * task.xdata.mean;
end

% Make Hessian exactly symmetric
q_posterior_cov_inv  = 0.5 * [q_posterior_cov_inv + q_posterior_cov_inv'];

q_posterior.mean     = q_posterior_cov_inv \ TT;


% ------------------------------------------------------------
% Posterior mode (under constraints)

xconstraints = [-task.xlower.value_nat; task.xupper.value_nat;];
Qconstraints = [-task.Q_xlower_q;   task.Q_xupper_q;];

if options.use_bounds_from_quantity_table,
  xconstraints = [xconstraints; -task.xall.lower_nat; task.xall.upper_nat;];
  Qconstraints = [Qconstraints; -task.Q_xall_q;       task.Q_xall_q];
end
  
q_posterior.mode = q_posterior.mean;

if flag_bounds, 
if sum([Qconstraints * q_posterior.mean > xconstraints]), 

  epsilon = 10^-10;

  %% Check whether constraints are already satisfied:
  %% [Qconstraints * q_posterior.mean < xconstraints - epsilon]

  lb = []; %lb = -10^15*ones(size(q_posterior.mean));
  ub = []; %ub =  10^15*ones(size(q_posterior.mean));

  if exist('cplexqp','file'),
    opt =  cplexoptimset('Display','off');
    [q_posterior.mode,fval,exitflag] = cplexqp(full(q_posterior_cov_inv), full(-q_posterior_cov_inv * q_posterior.mean), full(Qconstraints), xconstraints - epsilon,[],[],lb,ub,[],opt);
  else,
     opt = optimset('Display','off','Algorithm','active-set','MaxIter',10^8);
     [q_posterior.mode,fval,exitflag] = quadprog(full(q_posterior_cov_inv), full(-q_posterior_cov_inv * q_posterior.mean), full(Qconstraints), xconstraints - epsilon,[],[],lb,ub,[],opt);
  end

  if exitflag <0, error(sprintf('Error in optimisation during parameter balancing - exitflag %d',exitflag)); end

if exitflag~=1, 
  if exitflag == 0, 
    exitflag
    error('No solution found in quadratic programming; increase iteration number');
  end
  if exitflag == 5, 
    warning('Optimum found in flat region .. result may be unreliable');  
  else,
    exitflag
    warning('Problem in quadratic programming; possible reasons are too tight constraints (for concentrations, reaction affinities etc, or an infeasible flux distribution');  
  end
end

end
end

%  ----------------------------------------------------------

q_posterior.cov          = inv(q_posterior_cov_inv);
q_posterior.std          = sqrt(diag(q_posterior.cov));

xmodel_posterior.mode    = task.Q_xmodel_q * q_posterior.mode;
xmodel_posterior.mean    = task.Q_xmodel_q * q_posterior.mean;
xmodel_posterior.cov     = task.Q_xmodel_q * q_posterior.cov * task.Q_xmodel_q';

if isfield(options,'fix_Keq_in_sampling'),
if options.fix_Keq_in_sampling,
  display('  Keeping equilibrium constants fixed while computing the posterior covariance, standard deviations, and samples');
  my_ind = [task.q.indices.KV; task.q.indices.KM; task.q.indices.KA; task.q.indices.KI];
  xmodel_posterior.cov     = task.Q_xmodel_q(:,my_ind) * q_posterior.cov(my_ind,my_ind) * task.Q_xmodel_q(:,my_ind)';
end
end

xmodel_posterior.std     = sqrt(diag(full(xmodel_posterior.cov)));
q_posterior.samples      = repmat(q_posterior.mode,1,options.n_samples) ...
                           + real(sqrtm(full(q_posterior.cov))) * randn(length(q_posterior.mean),options.n_samples);
xmodel_posterior.samples = task.Q_xmodel_q * q_posterior.samples;

% [Qconstraints * q_posterior.samples, xconstraints - epsilon]
% prod(double(Qconstraints * q_posterior.samples <= repmat(xconstraints,1,options.n_samples)))

if length(task.xdata.mean),
  xdata_posterior.mode   = task.Q_xdata_q * q_posterior.mode;
  xdata_posterior.mean   = task.Q_xdata_q * q_posterior.mean;
  xdata_posterior.cov    = task.Q_xdata_q * q_posterior.cov * task.Q_xdata_q';
  xdata_posterior.std    = sqrt(diag(xdata_posterior.cov));
end

result.q_posterior.mode      = q_posterior.mode     ;
result.q_posterior.mean      = q_posterior.mean     ;
result.q_posterior.cov       = q_posterior.cov      ;
result.q_posterior.std       = q_posterior.std      ;
result.xmodel_posterior.mode = xmodel_posterior.mode;
result.xmodel_posterior.mean = xmodel_posterior.mean;
result.xmodel_posterior.cov  = xmodel_posterior.cov ;
result.xmodel_posterior.std  = xmodel_posterior.std ;
result.xmodel_posterior.samples = xmodel_posterior.samples;

if length(task.xdata.mean),
  result.xdata_posterior.mode = xdata_posterior.mode ;
  result.xdata_posterior.mean = xdata_posterior.mean ;
  result.xdata_posterior.cov  = xdata_posterior.cov  ;
  result.xdata_posterior.std  = xdata_posterior.std  ;
end

% -----------------------------------------------------------------------
% translate the posterior mean values back to kinetic constants


[nr,nm,nx,ind_KM,ind_KA,ind_KI,nKM,nKA,nKI] = network_numbers(task.network);

num_model = parameter_balancing_quantity_numbers(task.model_quantities,parameter_prior,task.network);

for it = 1:length(task.model_quantities),
  
  my_quantity     = task.model_quantities{it};
  ind             = find(strcmp(my_quantity,parameter_prior.QuantityType));
  my_scaling      = parameter_prior.MathematicalType{ind};
  my_symbol       = parameter_prior.Symbol{ind};
  my_indices      = task.xmodel.indices.(my_symbol);
  my_x_mode       = xmodel_posterior.mode(my_indices);
  my_x_mean       = xmodel_posterior.mean(my_indices);
  my_x_std        = xmodel_posterior.std(my_indices);
  my_x_samples    = xmodel_posterior.samples(my_indices,:);

  if strcmp(my_scaling, 'Multiplicative'),
     my_x_mode            = exp(my_x_mode);
     my_x_median          = exp(my_x_mean);
     my_x_geom_std        = exp(my_x_std);
     my_x_samples         = exp(my_x_samples);
    [my_x_mean, my_x_std] = lognormal_log2normal(my_x_mean, my_x_std,'arithmetic');
  else
    my_x_median          = nan * my_x_mode;
    my_x_geom_std        = nan * my_x_mode;
  end

 %figure(2); clf;
 %try
 %  plot(my_x_std,std(my_x_samples')','r.'); axis equal; pause
 %end
  
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
      dum_geom_std = dum; dum_geom_std(indices) = my_x_geom_std; my_x_geom_std = dum_geom_std;
      dum_samples = zeros(nr,nm,options.n_samples);
      for itt=1:options.n_samples,
	my_dum_samples = sparse(zeros(nr,nm));
	my_dum_samples(indices) = my_x_samples(:,itt);
	dum_samples(:,:,itt)   = my_dum_samples;
      end
      my_x_samples = dum_samples;
  end

  kpm_mode.(my_symbol)     = my_x_mode;
  kpm_mean.(my_symbol)     = my_x_mean;
  kpm_std.(my_symbol)      = my_x_std;
  kpm_median.(my_symbol)   = my_x_median;
  kpm_geom_std.(my_symbol) = my_x_geom_std;

  for itt=1:options.n_samples,
  switch my_symbol, 
    case {'KM', 'KA', 'KI'},
      kpm_samples{itt}.(my_symbol) = my_x_samples(:,:,itt);
    otherwise,
      kpm_samples{itt}.(my_symbol) = my_x_samples(:,itt);
    end
  end

end

% check thermodynamic constraints (in case of themrodynamics parameter balancing)
%
% [task.network.N' * kpm_mode.mu0 / RT, -log(kpm_mode.Keq)]
% [task.network.N' *[kpm_mode.mu0 + RT * log(kpm_mode.c)], -kpm_mode.A]


% -------------------------------------------------------------


result.kinetics.posterior_mode                   = kpm_mode;
result.kinetics.unconstrained_posterior_mean     = kpm_mean;
result.kinetics.unconstrained_posterior_std      = kpm_std;
result.kinetics.unconstrained_posterior_median   = kpm_median;
result.kinetics.unconstrained_posterior_geom_std = kpm_geom_std;

if options.n_samples,
  result.kinetics_posterior_samples = kpm_samples;
end

