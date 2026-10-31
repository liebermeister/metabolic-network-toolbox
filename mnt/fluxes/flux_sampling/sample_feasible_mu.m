function [mu, success_flag] = sample_feasible_mu(N,ind_ext,v,constraints,options,method,n_sample)

% SAMPLE_FEASIBLE_MU Compute chemical potentials that agree with a flux vector
%
% [mu, success_flag] = sample_feasible_mu(N,ind_ext,v,constraints,options,method,n_sample)
%
% Compute mu vectors that agree with a given flux vector v: delta_mu vectors in the feasible
% region (set by the previously sampled or determined flux directions and by bounds)
% are converted to mu vectors using the pseudoinverse of N'
%
% method: 'sample' or 'hit-and-run' (default): sample n_sample dmu vectors uniformly from the
%                           feasible polytope (by hit-and-run sampling)
% method: 'extreme_points': vertices of the feasible polytope with minimal and maximal values
%                           along the coordinate axes (this is no actual sampling)
% method: 'centre' or 'central': central dmu vector
%
% See also: sample_feasible_v, sample_feasible_dmu

eval(default('method','''sample''','n_sample','1','options','struct'));

mu           = [];
success_flag = 1;

[nm,nr] = size(N);

if ~isfield(options,'verbose'),         options.verbose = 0; end
if ~isfield(options,'seed'),            options.seed = 0; end
if ~isnan(options.seed), rng(options.seed); end

% compute minimal bounding box for reduced chemical potentials
% U consists of the linearly independent columns of N'
% represent delta_mu as U * mured
% represent       N = D * U'

[ech,ind] = rref(N);
U         = ech(1:length(ind),:)';
n_mured   = size(U,2);
D         = round(N * pinv(full(U')));

% WORKAROUND
constraints.mured_min = mean(constraints.mu_min)*ones(n_mured,1);
constraints.mured_max = mean(constraints.mu_max)*ones(n_mured,1);

dmu_limit      = constraints.dmu_limit - 10^-5; % be a little stricter
dmu_limit_min  = constraints.dmu_limit_min + 10^-5; % be a little stricter

constraints.dmu_min(constraints.dmu_min<-dmu_limit) = -dmu_limit;
constraints.dmu_max(constraints.dmu_max> dmu_limit) =  dmu_limit;

constraints.dmu_max(find(double([v>0]).*double([constraints.dmu_max>-dmu_limit_min]))) = -dmu_limit_min;
constraints.dmu_min(find(double([v<0]).*double([constraints.dmu_min< dmu_limit_min]))) =  dmu_limit_min;

ind_dmu_fix = find(isfinite(constraints.dmu_fix));

if ind_dmu_fix
  display('Using fixed dmu values'),
  constraints.dmu_max(ind_dmu_fix) = constraints.dmu_fix(ind_dmu_fix);
  constraints.dmu_min(ind_dmu_fix) = constraints.dmu_fix(ind_dmu_fix);
end

epsilon     = RT/constraints.rho;
v_signs     = sign(v);
ind_nonzero = find(v_signs~=0);
v_sign_U    = diag(v_signs) * U;
v_sign_U    = v_sign_U(ind_nonzero,:);

A =  [- eye(n_mured); eye(n_mured); ...
     [- eye(nr); eye(nr) ] * U;
        v_sign_U];

b = [- constraints.mured_min; constraints.mured_max;
     - constraints.dmu_min;   constraints.dmu_max;...
     - epsilon*ones(length(ind_nonzero),1)];

% avoid inf values in linprog
b(isinf(b))= sign(b(isinf(b))) * 10^10;

switch method

  case 'extreme_points'

    %% optimise: vec * mured = min while A * mured <= b
    %%
    %% that is, in detail:
    %% - mured             <= mured_min   i.e.  mured >= mured_min
    %%   mured             <= mured_max   i.e.  mured <= mured_miax
    %% - U mured           <= dmu_min     i.e.  dmu   >= dmu_min
    %%   U mured           <= dmu_max     i.e.  dmu   <= dmu_max
    %%  dg(v_sign)*U*mured <= -epsilon    i.e.  sign(v) * dmu < epsilon <= 0

    for it = 1:n_mured
      if options.verbose, display(sprintf('Dimension %d/%d',it,n_mured)); end
      vec = zeros(n_mured,1); vec(it) = 1;
      [mured_lower_opt,fval,exitflag] = lp_solve(vec, A, b);
      if exitflag <=0, exitflag
        error('Error during linear programming problem for finding feasible chemical potentials');
      end
      [mured_upper_opt,fval,exitflag] = lp_solve(-vec, A, b);
      if exitflag <=0, exitflag
        error('Error during linear programming problem for finding feasible chemical potentials');
      end
      dmu_lower(:,it) = U * mured_lower_opt;
      dmu_upper(:,it) = U * mured_upper_opt;
    end

    all_dmu = [dmu_lower, dmu_upper];

  case {'centre','central'}
    all_dmu = U * find_polytope_centre([],[],A,b,[],[]);

  case {'sample','hit-and-run'}
    [all_dmu, success_flag] = sample_dmu_hit_and_run(U,A,b,ind_dmu_fix,constraints.dmu_fix,n_sample);
    if ~success_flag, return; end

  otherwise
    error('Unknown method %s',method);

end

mu = pinv(full(N')) * all_dmu;

% ------------------------------------------------------------------------------

function [dmu, success_flag] = sample_dmu_hit_and_run(U,A,b,ind_fix,dmu_fix,n_sample)

% Sample mured uniformly from the polytope A * mured <= b. Since U has full column rank,
% dmu = U * mured is uniform as well. Predefined dmu values make the polytope flat, so
% they are treated as equality constraints: mured = mured_fix + K * y

success_flag = 1;
dmu          = [];
n_mured      = size(U,2);
tolerance    = 10^-8;

if isempty(ind_fix)
  mured_fix = zeros(n_mured,1);
  K         = eye(n_mured);
else
  mured_fix = pinv(U(ind_fix,:)) * dmu_fix(ind_fix);
  K         = null(U(ind_fix,:));
  if max(abs(U(ind_fix,:) * mured_fix - dmu_fix(ind_fix))) > tolerance
    success_flag = 0;
    return
  end
end

M  = A * K;
bb = b - A * mured_fix;

% omit constraints that do not depend on y, and the formal bounds 10^10 (see above)
is_constant = sqrt(sum(M.^2,2)) < tolerance;
if sum(bb(is_constant) < -tolerance)
  success_flag = 0;
  return
end
keep = find(~is_constant .* (bb < 10^9));
M    = M(keep,:);
bb   = bb(keep);

n_y = size(K,2);

if n_y == 0
  dmu = U * repmat(mured_fix,1,n_sample);
  return
end

y_start = find_polytope_centre([],[],M,bb,[],[]);
n_warm  = 1000 + 100 * n_y;
n_pick  = 10 + n_y;
y       = convex_sampling_hit_and_run(M,bb,[],[],[],[],n_sample,n_warm,n_pick,0,y_start,0);
dmu     = U * (repmat(mured_fix,1,n_sample) + K * y);
