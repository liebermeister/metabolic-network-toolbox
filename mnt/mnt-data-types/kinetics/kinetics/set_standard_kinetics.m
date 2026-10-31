% kinetics = set_standard_kinetics(network,parameters)
%
% construct a kinetics field for a metabolic network
% set default kinetics types for all reactions of a network,
% depending on number of substrates and products and reversibility
%
% kinetics         parameters
%
% 'rev_uni_uni'    {'k','km','Pk','Pkm'}
% 'rev_uni_bi'     {'k','km','Pk','Pkm'}
% 'rev_bi_uni'     {'k','km','Pk','Pkm'}
% 'rev_bi_bi'      {'k','km','Pk','Pkm'}
% 'irrev_uni'      {'k','Pk'}
% 'irrev_bi'       {'k','Pk'}
% 'influx'         {'v_in'}
% 'efflux'         {'k_out'}
% 'rev_multiple'   {'k','km','Pk','Pkm'}

function kinetics = set_standard_kinetics(network,parameters)

if ~exist('parameters','var'), parameters = struct('ready_made',1); end

kinetics.type      = 'standard';
kinetics.reactions = cell(length(network.actions),1);

for i=1:length(network.actions)

  kinetics.reactions{i}.type = 'unknown';
  rev   = network.reversible(i);
  n_in  = sum(network.N(:,i)<0) ;
  n_out = sum(network.N(:,i)>0) ;
  n_act = 0;
  n_rep = 0;
  if isfield(network,'regulation_matrix'),
    n_act = full(sum(network.regulation_matrix(i,:)>0));
    n_rep = full(sum(network.regulation_matrix(i,:)<0));
  end

%  switch rev
%    case 1,
%      if n_in == 1 & n_out == 1, kinetics.reactions{i}.type = 'rev_uni_uni'; end
%      if n_in == 1 & n_out == 2, kinetics.reactions{i}.type = 'rev_uni_bi'; end
%      if n_in == 2 & n_out == 1, kinetics.reactions{i}.type = 'rev_bi_uni'; end
%      if n_in == 2 & n_out == 2, kinetics.reactions{i}.type = 'rev_bi_bi'; end
%    case 0,
%      if n_in == 1 , kinetics.reactions{i}.type = 'irrev_uni'; end
%      if n_in == 2 , kinetics.reactions{i}.type = 'irrev_bi'; end
%  end

  switch rev
    case 1,
      if n_in == 1 & n_out == 1, kinetics.reactions{i}.type = 'standard_rev_uni_uni'; end
      if n_in == 1 & n_out == 2, kinetics.reactions{i}.type = 'standard_rev_uni_bi'; end
      if n_in == 2 & n_out == 1, kinetics.reactions{i}.type = 'standard_rev_bi_uni'; end
      if n_in == 2 & n_out == 2, kinetics.reactions{i}.type = 'standard_rev_bi_bi'; end
    case 0,
      if n_in == 1 , kinetics.reactions{i}.type = 'standard_irrev_uni'; end
      if n_in == 2 , kinetics.reactions{i}.type = 'standard_irrev_bi'; end
  end

%  switch rev,
%    case 1,
%      kinetics.reactions{i}.type='multiple_michaelis_menten_rev';
%    case 0,
%      kinetics.reactions{i}.type='multiple_michaelis_menten_irrev';
%  end

  if n_in == 0 & n_out > 0,  kinetics.reactions{i}.type = 'influx'; end
  if n_in > 0  & n_out == 0, kinetics.reactions{i}.type = 'efflux'; end
  if strcmp(kinetics.reactions{i}.type,'unknown'),
  switch rev
    case 1,
      kinetics.reactions{i}.type='standard_rev_multiple';
    case 0,
      kinetics.reactions{i}.type='irrev_multiple';
  end
  end

  if isfield(parameters,'ready_made'), if parameters.ready_made==1, kinetics.reactions{i}.type='ready_made'; end; end

  switch kinetics.reactions{i}.type,

    case 'influx',
      kinetics.reactions{i}.parameters = {'v_in'}';
      kinetics.reactions{i}.sizes = [1]';
    case 'efflux',
      kinetics.reactions{i}.parameters = {'k_out'}';
      kinetics.reactions{i}.sizes = [1]';

    case 'mass-action',
      kinetics.reactions{i}.parameters = {'k','km'}';
      kinetics.reactions{i}.sizes = [1 1]';

    case 'michaelis_menten_rev',
      kinetics.reactions{i}.parameters = {'Km_S','Km_P','Vm_S','Vm_P'}';
      kinetics.reactions{i}.sizes = [1 1 1 1]';
    case 'multiple_michaelis_menten_rev',
      kinetics.reactions{i}.parameters = {'Km_S','Km_P','Vm_S','Vm_P'}';
      kinetics.reactions{i}.sizes = [1 1 1 1]';
    case 'multiple_michaelis_menten_irrev',
      kinetics.reactions{i}.parameters = {'Km_S','Vm_S'}';
      kinetics.reactions{i}.sizes = [1 1]';
    case {'irrev_uni'},
      kinetics.reactions{i}.parameters = {'k','Pk'}';
      kinetics.reactions{i}.sizes = [2 1]';
    case {'irrev_bi'},
      kinetics.reactions{i}.parameters = {'k','Pk'}';
      kinetics.reactions{i}.sizes = [3 1]';
    case {'irrev_multiple'},
      kinetics.reactions{i}.parameters = {'k','Pk'}';
      kinetics.reactions{i}.sizes = [2 1]';
    case {'rev_uni_uni'},
      kinetics.reactions{i}.parameters = {'k','km','Pk','Pkm'}';
      kinetics.reactions{i}.sizes = [2 2 1 1]';
    case {'rev_uni_bi'},
      kinetics.reactions{i}.parameters = {'k','km','Pk','Pkm'}';
      kinetics.reactions{i}.sizes = [3 3 1 1]';
    case {'rev_bi_uni'},
      kinetics.reactions{i}.parameters = {'k','km','Pk','Pkm'}';
      kinetics.reactions{i}.sizes = [3 3 1 1]';
    case {'rev_bi_bi'},
      kinetics.reactions{i}.parameters = {'k','km','Pk','Pkm'}';
      kinetics.reactions{i}.sizes = [4 4 1 1]';
    case {'rev_multiple'},
      kinetics.reactions{i}.parameters = {'k','km','Pk','Pkm'}';
      kinetics.reactions{i}.sizes = [2 2 1 1]';

    case {'standard_irrev_uni'},
      kinetics.reactions{i}.parameters = {'Vf','kA'}';
      kinetics.reactions{i}.sizes = [1 1]';
    case {'standard_irrev_bi'},
      kinetics.reactions{i}.parameters = {'Vf','kA','kB'}';
      kinetics.reactions{i}.sizes = [1 1 1]';
    case {'standard_rev_uni_uni'},
      kinetics.reactions{i}.parameters = {'Vf','Vr','kA','kP'}';
      kinetics.reactions{i}.sizes = [1 1 1 1]';
    case {'standard_rev_uni_bi'},
      kinetics.reactions{i}.parameters = {'Vf','Vr','kA','kP','kQ'}';
      kinetics.reactions{i}.sizes = [1 1 1 1 1]';
    case {'standard_rev_bi_uni'},
      kinetics.reactions{i}.parameters = {'Vf','Vr','kA','kB','kP'}';
      kinetics.reactions{i}.sizes = [1 1 1 1 1]';
    case {'standard_rev_bi_bi'},
      kinetics.reactions{i}.parameters = {'Vf','Vr','kA','kB','kP','kQ'}';
      kinetics.reactions{i}.sizes = [1 1 1 1 1 1]';
    case {'standard_rev_multiple'},
      kinetics.reactions{i}.parameters = {'Vf','Vr','kf','kr'}';
      kinetics.reactions{i}.sizes = full([1 1 n_in n_out])';

    case {'ready_made'},
      kinetics.reactions{i}.parameters = {'Vf','Vr','KMf','KMr','KA','KI'}';
      kinetics.reactions{i}.sizes = full([1 1 n_in n_out n_act n_rep])';

  end

  switch kinetics.reactions{i}.type,
    case {'standard_irrev_uni','standard_irrev_bi','standard_rev_uni_uni','standard_rev_uni_bi','standard_rev_bi_uni','standard_rev_bi_bi','standard_rev_multiple'},
     if n_act > 0,
       kinetics.reactions{i}.parameters = [kinetics.reactions{i}.parameters; {'Kact'} ];
       kinetics.reactions{i}.sizes = [kinetics.reactions{i}.sizes; n_act];
     end
     if n_rep > 0,
       kinetics.reactions{i}.parameters = [kinetics.reactions{i}.parameters; {'Krep'} ];
       kinetics.reactions{i}.sizes = [kinetics.reactions{i}.sizes; n_rep];
     end
  end

  for it=1:length(kinetics.reactions{i}.parameters),
    kinetics.reactions{i}=setfield(kinetics.reactions{i},kinetics.reactions{i}.parameters{it}, ...
                                                 ones(kinetics.reactions{i}.sizes(it),1));
  end
end
