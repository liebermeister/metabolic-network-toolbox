% kin = vectorchange2standard_kinetics(kinetics,indices_changed,offsets,reaction_pars,par_per)
%
%Construct kinetics structure from a parameter vector (for kinetics of type 'standard')

function [kinetics,it] = vectorchange2standard_kinetics(kinetics,indices_changed,offsets,reaction_pars,par_per)

it = sum( indices_changed > offsets );

% reaction number

parnum    = indices_changed-offsets(it);
v         = reaction_pars{it};
v(parnum) = par_per;
kinetics.reactions{it} = set_reaction(kinetics.reactions{it},v);

% -------------------------------------------------------------------------------

function rea = set_reaction(rea,v);

switch rea.type,

  case 'ready_made',
    rea.Vf  = v(1);
    rea.Vr  = v(2);
    rea.KMf = v(2 + (1:rea.sizes(3)));
    rea.KMr = v(2 + rea.sizes(3) + (1:rea.sizes(4)));
    rea.KA  = v(2 + sum(rea.sizes(3:4)) + (1:rea.sizes(5)));
    rea.KI  = v(2 + sum(rea.sizes(3:5)) + (1:rea.sizes(6)));

  otherwise,

    switch rea.type,
      case 'mass-action',                      rea.k=v(1); rea.km=v(2);
      case 'mass-action_rev',                  rea.k=v(1); rea.km=v(2);
      case 'mass-action_irrev',                rea.k=v(1);
      case 'multiple_michaelis_menten_rev',    rea.Km_S=v(1); rea.Km_P=v(2); rea.Vm_S=v(3); rea.Vm_P=v(4);
      case 'multiple_michaelis_menten_irrev',  rea.Km_S=v(1); rea.Vm_S=v(2);
      case 'michaelis_menten_rev',             rea.Vm_S=v(1);  rea.Km_S=v(2); rea.Vm_P=v(3);  rea.Km_P=v(4);
      case 'influx',                rea.v_in  = v(1);
      case 'efflux',                rea.k_out = v(1);
      case 'rev_uni_uni',           rea.Pk=v(1);  rea.Pkm=v(2); rea.k=v(3);  rea.km=v(4);
      case 'rev_bi_uni',            rea.Pk=v(1);  rea.Pkm=v(2); rea.k=v(3);  rea.km=v(4);
      case 'rev_uni_bi',            rea.Pk=v(1);  rea.Pkm=v(2); rea.k=v(3);  rea.km=v(4);
      case 'rev_bi_bi',             rea.Pk=v(1);  rea.Pkm=v(2); rea.k=v(3);  rea.km=v(4);
      case 'irrev_uni',             rea.Pk=v(1);  rea.k=v(2);
      case 'irrev_bi',              rea.Pk=v(1);  rea.k=v(2);
      case 'irrev_multiple',        rea.k  = v(1:2);  rea.Pk = v(3);
      case 'rev_multiple',          rea.k  = v(1:2); rea.km = v(3:4); rea.Pk = v(5); rea.Pkm = v(6);
      case 'standard_irrev_uni',    rea.Vf = v(1); rea.kA=v(2);
      case 'standard_irrev_bi',     rea.Vf = v(1); rea.kA=v(2); rea.kB=v(3);
      case 'standard_rev_uni_uni',  rea.Vf = v(1); rea.Vr=v(2); rea.kA=v(3); rea.kP=v(4);
      case 'standard_rev_uni_bi',   rea.Vf = v(1); rea.Vr=v(2); rea.kA=v(3); rea.kP=v(4); rea.kQ=v(5);
      case 'standard_rev_bi_uni',   rea.Vf = v(1); rea.Vr=v(2); rea.kA=v(3); rea.kB=v(4); rea.kP=v(5);
      case 'standard_rev_bi_bi',    rea.Vf = v(1); rea.Vr=v(2); rea.kA=v(3); rea.kB=v(4); rea.kP=v(5); rea.kQ=v(6);
      case 'standard_rev_multiple', rea.Vf = v(1); rea.Vr=v(2); rea.kf=v(3:2+rea.sizes(3)); ...
            rea.kr=v(3+rea.sizes(3):2+rea.sizes(3)+rea.sizes(4));
      otherwise error('unknown kinetics type');
    end

    nrep = 0;
    nact = 0;
    if isfield(rea,'Krep'),   nrep = length(rea.Krep);    end
    if isfield(rea,'Kact'),   nact = length(rea.Kact);    end

    if isfield(rea,'Krep'),
      rea.Krep = v(sum(rea.sizes)+1-nrep-nact:sum(rea.sizes)-nact);
    end
    if isfield(rea,'Kact'),
      rea.Kact = v(sum(rea.sizes)+1-nact:sum(rea.sizes));
    end
end