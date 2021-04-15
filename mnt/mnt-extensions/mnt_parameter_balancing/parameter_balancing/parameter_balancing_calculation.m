function [result, exitflag] = parameter_balancing_calculation(task, parameter_prior, pb_options)
  
% [result, exitflag] = parameter_balancing(task, parameter_prior, pb_options)
%
% Parameter balancing based on vectors and matrices in data structure 'task' 
%
% Function inputs (mandatory)
%  task:  parameter balancing task (see parameter_balancing_task.m)
%
% Function inputs (optional)
%  parameter_prior: table of biochemical quantities (see biochemical_parameter_prior.m)
%  pb_options:       options table with fields 
%                   pb_options.use_pseudo_values (default 0)
%                   pb_options.n_samples            (default 0)
%                     number of random samples from posterior 
%                     for output result.kinetics_posterior_samples
%                   pb_options.fix_Keq_in_sampling

global log_text % text for the log file is added to this variable 

eval(default('parameter_prior','[]','pb_options','struct'));

if isempty(parameter_prior),
  parameter_prior = biochemical_parameter_prior;
end

pb_options = join_struct(parameter_balancing_options,pb_options);

q_prior.cov_inv = diag(sparse(1./task.q.prior.std.^2));

if find(task.xdata.std==0),
  warning('Vanishing standard deviations found. Replacing them by 10^-3.'); 
  task.xdata.std(task.xdata.std==0) = 10^-3;
end

if find(task.xdata.std<10^-3),
  warning('Vanishing standard deviations found. Replacing them by 10^-3.'); 
  task.xdata.std(task.xdata.std<10^-3) = 10^-3;
end

xdata.cov_inv = diag(sparse(1./task.xdata.std.^2));

% --------------------
% Use predefined covariance matrix for Delta G of reaction (if provided in file pb_options.dg_cov_file)
% --------------------

if isfield(task.xdata.indices,'dmu0') * isfield(pb_options,'dg_precision_file'),
  if length(pb_options.dg_precision_file),
    display(sprintf('Reading Delta G covariance matrix from file %s',pb_options.dg_precision_file));
    load(pb_options.dg_precision_file);
    % if rxn_id is not a cell array (but a character array), then convert it to a cell array
    if ischar(rxn_id),
      rxn_id = cellstr(rxn_id);
      ind_reactions_in_dg_precision_matrix = label_names(rxn_id,task.network.actions);
      % Now there are a few tests to make sure that the covariance matrix agrees with the dmu0 data 
      % indices must be increasing
      if norm(ind_reactions_in_dg_precision_matrix-sort(ind_reactions_in_dg_precision_matrix)),
        error('data in Delta G covariance matrix must have the same order as in the model');
      end
      % size of covariance matrix must match the number of dmu0 data points
      if length(ind_reactions_in_dg_precision_matrix) ~= length(task.xdata.std(task.xdata.indices.dmu0)),
        error('Size of Delta G covariance matrix must match the number of data points. Maybe some data points have been filtered out because they were outside the allowed bounds? This must be fixed.');
      end
      %% diagonal values of covariance must match the standard deviations in the data
      %% This test ist discarded because with poorly defined
      %% precisions in some directions, the comparison may be erroneous
      % if max(abs([task.xdata.std(task.xdata.indices.dmu0)./sqrt(diag(dg_precision))]-1))>0.02,
      %   [task.xdata.std(task.xdata.indices.dmu0)./sqrt(diag(dg_precision))]
      %   error('The diagonal values of the Delta G covariance matrix must agree with the standard deviations in the Delta G data.');
      % end
      % check the eigenvalues
      eig_min = min(eig(dg_precision));
      eig_max = max(eig(dg_precision));
      if eig_min < (10^-10)*eig_max,
        display('  (parameter_balancing_calculation.m): The delta G precision matrix contains negative eigenvalues or is badly conditioned. I shift all eigenvalues up');
        dg_precision = dg_precision + [-min(eig_min,0)+(10^-10)*eig_max] * eye(size(dg_precision,1));
      end
      display('Using prepared covariance matrix for standard delta G values');
      xdata.cov_inv(task.xdata.indices.dmu0,task.xdata.indices.dmu0) = dg_precision;
    end
  end
end

xall_pseudo.cov_inv = diag(sparse(1./task.xpseudo.std.^2));


% ------------------------------------------------------------
% Compute posterior mean and precision matrix (without constraints)
% q_posterior.mean and q_posterior_cov_inv
% ------------------------------------------------------------

% Start with prior terms
q_posterior_cov_inv = q_prior.cov_inv;
TT                  = q_prior.cov_inv * task.q.prior.mean;

%% Add pseudo value terms
if pb_options.use_pseudo_values,
  if prod(size(xall_pseudo.cov_inv)),
    q_posterior_cov_inv = q_posterior_cov_inv + task.Q_xpseudo_q' * xall_pseudo.cov_inv * task.Q_xpseudo_q;
    TT                  = TT                  + task.Q_xpseudo_q' * xall_pseudo.cov_inv * task.xpseudo.mean;
  end
end

% Add data terms
if length(task.Q_xdata_q),
  q_posterior_cov_inv = q_posterior_cov_inv + task.Q_xdata_q' * xdata.cov_inv * task.Q_xdata_q;
  TT                  = TT                  + task.Q_xdata_q' * xdata.cov_inv * task.xdata.mean;
end

% Make Hessian exactly symmetric
q_posterior_cov_inv  = 0.5 * [q_posterior_cov_inv + q_posterior_cov_inv'];

% Check for ill-conditioned inverse posterior covariance matrix
eigenvalues = sort(eig(full(q_posterior_cov_inv))); 

if find(eigenvalues<0), 
  error('Inverse posterior covariance matrix has negative eigenvalues - please check for unrealistically small error bars in your data file');
end

if abs(max(eigenvalues) / min(eigenvalues)) > 10^15,
  log_text = [log_text, '\n  (parameter_balancing_calculation.m): Inverse posterior covariance matrix appears ill-conditioned - please check for unrealistically small error bars in your data file'];
  max_eigenvalue = max(eigenvalues)
  min_eigenvalue = min(eigenvalues)
  figure(1000); plot(eigenvalues); 
  xlabel('ordering'); ylabel('eigenvalue of inverse posterior covariance matrix'); set(gca,'YScale','Log');
end

q_posterior.mean = q_posterior_cov_inv \ TT;


% ------------------------------------------------------------
% Compute posterior mode q_posterior.mode (under constraints)
% ------------------------------------------------------------

% --- build vectors and matrices for equality and inequality constraints

ind_eq = find(task.xlower.value_nat == task.xupper.value_nat);
ind_in = find(task.xlower.value_nat < task.xupper.value_nat);
if find(task.xlower.value_nat > task.xupper.value_nat), error('Infeasible bounds'); end

eqconstraints_names = task.xlower.names(ind_eq);
xeqconstraints      = task.xlower.value_nat(ind_eq);
Qeqconstraints      = task.Q_xlower_q(ind_eq,:);

inconstraints_names = [ task.xlower.names(ind_in);     task.xupper.names(ind_in)];
xconstraints        = [-task.xlower.value_nat(ind_in); task.xupper.value_nat(ind_in)];
Qconstraints        = [-task.Q_xlower_q(ind_in,:);     task.Q_xupper_q(ind_in,:)];

% --- use the constraints?

use_constraints = length(task.xlower.value_nat) + length(task.xupper.value_nat) > 0;

if pb_options.ignore_all_constraints,
  log_text = [log_text, '\no Ignoring all constraints (option "ignore_all_constraints" was chosen)'];
  use_constraints = 0;
end

if use_constraints,
  % make the constraints a bit stricter, to ensure that solutions (even with numerical errors) will satisfy the constraints 
  epsilon = 10^-4;
  %% Project solution to feasible point (satisfying constraints) (using function at the bottom of this file)
  q_posterior.mode = compute_posterior_mode_with_constraints(q_posterior.mean, q_posterior_cov_inv, Qconstraints, xconstraints, Qeqconstraints, xeqconstraints, epsilon);
else
  q_posterior.mode = q_posterior.mean;
end

%  ----------------------------------------------------------

q_posterior.cov          = inv(q_posterior_cov_inv);
q_posterior.std          = sqrt(diag(q_posterior.cov));

xmodel_posterior.mode    = task.Q_xmodel_q * q_posterior.mode;
xmodel_posterior.mean    = task.Q_xmodel_q * q_posterior.mean;
xmodel_posterior.cov     = task.Q_xmodel_q * q_posterior.cov * task.Q_xmodel_q';

if isfield(pb_options,'fix_Keq_in_sampling'),
  if pb_options.fix_Keq_in_sampling,
    log_text = [log_text, '\n  Keeping equilibrium constants fixed while computing the posterior covariance, standard deviations, and samples'];
    my_ind = [task.q.indices.KV; task.q.indices.KM; task.q.indices.KA; task.q.indices.KI];
    xmodel_posterior.cov     = task.Q_xmodel_q(:,my_ind) * q_posterior.cov(my_ind,my_ind) * task.Q_xmodel_q(:,my_ind)';
  end
end

xmodel_posterior.std     = sqrt(diag(full(xmodel_posterior.cov)));

if pb_options.n_samples >0,
  display(sprintf('o Generating %d samples from the posterior distribution', pb_options.n_samples ));
  q_posterior.samples = repmat(q_posterior.mode,1,pb_options.n_samples) + ...
      real(sqrtm(full(q_posterior.cov))) * randn(length(q_posterior.mean),pb_options.n_samples);
  
  if use_constraints,
    %% Project samples to feasible points (satisfying constraints) (using function at the bottom of this file)    
    for its = 1:pb_options.n_samples,
      q_posterior.samples(:,its) = compute_posterior_mode_with_constraints(q_posterior.samples(:,its), q_posterior_cov_inv, Qconstraints, xconstraints, Qeqconstraints, xeqconstraints, epsilon);
    end
    xmodel_posterior.samples = task.Q_xmodel_q * q_posterior.samples;
  end

else
  xmodel_posterior.samples = xmodel_posterior.mode;
end

% [Qconstraints * q_posterior.samples, xconstraints - epsilon]
% prod(double(Qconstraints * q_posterior.samples <= repmat(xconstraints,1,pb_options.n_samples)))

if length(task.xdata.mean),
  xdata_posterior.mode   = task.Q_xdata_q * q_posterior.mode;
  xdata_posterior.mean   = task.Q_xdata_q * q_posterior.mean;
  xdata_posterior.cov    = task.Q_xdata_q * q_posterior.cov * task.Q_xdata_q';
  xdata_posterior.std    = sqrt(diag(xdata_posterior.cov));
end


% -----------------------------------------------------------------------
% data structure 'result'

result.q_posterior.mode         = q_posterior.mode     ;
result.q_posterior.mean         = q_posterior.mean     ;
result.q_posterior.cov          = q_posterior.cov      ;
result.q_posterior.std          = q_posterior.std      ;
result.xmodel_posterior.mode    = xmodel_posterior.mode;
result.xmodel_posterior.mean    = xmodel_posterior.mean;
result.xmodel_posterior.cov     = xmodel_posterior.cov ;
result.xmodel_posterior.std     = xmodel_posterior.std ;
result.xmodel_posterior.samples = xmodel_posterior.samples;

result.constraints_on_q_Aineq  = Qconstraints; 
result.constraints_on_q_bineq  = xconstraints; 

if length(task.xdata.mean),
  result.xdata_posterior.mode = xdata_posterior.mode ;
  result.xdata_posterior.mean = xdata_posterior.mean ;
  result.xdata_posterior.cov  = xdata_posterior.cov  ;
  result.xdata_posterior.std  = xdata_posterior.std  ;
end


% -----------------------------------------------------------------------
% convert posterior mean values back to kinetic constants

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
      dum_samples = zeros(nr,nm,pb_options.n_samples);
      for itt=1:pb_options.n_samples,
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

  for itt=1:pb_options.n_samples,
  switch my_symbol, 
    case {'KM', 'KA', 'KI'},
      kpm_samples{itt}.(my_symbol) = myround(my_x_samples(:,:,itt));
    otherwise,
      kpm_samples{itt}.(my_symbol) = myround(my_x_samples(:,itt));
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

if pb_options.n_samples,
  result.kinetics_posterior_samples = kpm_samples;
end

% -------------------------------------------------------------------------------------


% -------------------------------------------------------------------------------------
% function for computing the posterior mode under constraints
% -------------------------------------------------------------------------------------

function q_posterior_mode = compute_posterior_mode_with_constraints(q_posterior_mean, q_posterior_cov_inv, Qconstraints, xconstraints, Qeqconstraints, xeqconstraints, epsilon)

global log_text % text for the log file is added to this variable 

% check whether q_posterior_mean violates the constraints

ind_active_constraints = find(Qconstraints * q_posterior_mean > xconstraints - epsilon);

if length(ind_active_constraints) == 0,
  q_posterior_mode = q_posterior_mean;

else
  lb = []; %lb = -10^15*ones(size(q_posterior_mean));
  ub = []; %ub =  10^15*ones(size(q_posterior_mean));

  if exist('cplexqp','file') * [size(q_posterior_cov_inv,1)>1000],
    display('  (parameter_balancing_calculation.m): Problem size is too large for CPLEX community edition; using matlab quadprog function instead')
  end

  if exist('cplexqp','file')* [size(q_posterior_cov_inv,1)<=1000],
    log_text = [log_text, 'Using CPLEX for quadratic optimisation'];
    opt =  cplexoptimset('Display','off');
    [q_posterior_mode,fval,exitflag] = cplexqp(full(q_posterior_cov_inv), full(-q_posterior_cov_inv * q_posterior_mean), full(Qconstraints), xconstraints - epsilon,full(Qeqconstraints), xeqconstraints,lb,ub,[],opt);
  else,
    log_text = [log_text, 'Using Matlab quadprog for quadratic optimisation'];
     opt = optimset('Display','off','Algorithm','interior-point-convex','MaxIter',10^8); % 'active-set'
     [q_posterior_mode,fval,exitflag] = quadprog(full(q_posterior_cov_inv), full(-q_posterior_cov_inv * q_posterior_mean), full(Qconstraints), xconstraints - epsilon,full(Qeqconstraints), xeqconstraints,lb,ub,[],opt);
  end

  if exitflag~=1,
    if exitflag < 0,  error(sprintf('Error in optimisation during parameter balancing - exitflag %d',exitflag)); end
    if exitflag == 0, error(sprintf('No solution found in quadratic programming (exitflag %d); increase iteration number',exitflag)); end
    if exitflag == 5, warning('Objective function around the optimum is flat; numerical result may be unreliable');  
    else,
      warning(sprintf('Problem in quadratic programming (existflag %d); maybe constraints are too tight (for concentrations, reaction affinities etc, or flux distribution is infeasible', exitflag));  
    end
  end

  % check 
  % check_ineq = prod([Qconstraints * q_posterior_mode <= xconstraints]);
  % check_eq   = prod([Qeqconstraints * q_posterior_mode == xeqconstraints]);
  
end

% -------------------------------------------------

function myround = myround(x)

myround = round(x,4,'significant');
