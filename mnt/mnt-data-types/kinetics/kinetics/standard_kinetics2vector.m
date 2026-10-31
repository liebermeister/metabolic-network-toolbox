%[v,names,relevant_reactions] = standard_kinetics2vector(kinetics)
%
%extract parameter vector from kinetics of type 'standard'

function [v,names,relevant_reactions] = standard_kinetics2vector(kinetics)
v=[];
names={};
if isfield(kinetics,'reactions'),
  r = kinetics.reactions;
  for it=1:length(r),
    switch r{it}.type,
      case 'mass-action',                      this_v = [ r{it}.k;    r{it}.km];
      case 'mass-action_rev',                  this_v = [ r{it}.k;    r{it}.km];
      case 'mass-action_irrev',                this_v = [ r{it}.k; ];
      case 'multiple_michaelis_menten_rev',    this_v = [ r{it}.Km_S; r{it}.Km_P;  r{it}.Vm_S;  r{it}.Vm_P];
      case 'multiple_michaelis_menten_irrev',  this_v = [ r{it}.Km_S; r{it}.Vm_S  ];
      case 'michaelis_menten_rev',             this_v = [ r{it}.Vm_S; r{it}.Km_S; r{it}.Vm_P;  r{it}.Km_P ];
      case 'influx',                           this_v = [  r{it}.v_in ];
      case 'efflux',                           this_v = [  r{it}.k_out];
      case 'rev_uni_uni',          this_v = [ r{it}.Pk;  r{it}.Pkm; r{it}.k;  r{it}.km];
      case 'rev_bi_uni',           this_v = [ r{it}.Pk;  r{it}.Pkm; r{it}.k;  r{it}.km];
      case 'rev_uni_bi',           this_v = [ r{it}.Pk;  r{it}.Pkm; r{it}.k;  r{it}.km];
      case 'rev_bi_bi',            this_v = [ r{it}.Pk;  r{it}.Pkm; r{it}.k;  r{it}.km];
      case 'irrev_uni',            this_v = [ r{it}.Pk;  r{it}.k];
      case 'irrev_bi',             this_v = [ r{it}.Pk;  r{it}.k];
      case 'irrev_multiple',       this_v = [ r{it}.k;   r{it}.Pk];
      case 'rev_multiple',         this_v = [ r{it}.k ;  r{it}.km; r{it}.Pk; r{it}.Pkm];
      case 'standard_irrev_uni',   this_v = [ r{it}.Vf;  r{it}.kA];
      case 'standard_irrev_bi',    this_v = [ r{it}.Vf;  r{it}.kA; r{it}.kB];
      case 'standard_rev_uni_uni', this_v = [ r{it}.Vf;  r{it}.Vr;  r{it}.kA; r{it}.kP];
      case 'standard_rev_uni_bi',  this_v = [ r{it}.Vf;  r{it}.Vr;  r{it}.kA; r{it}.kP; r{it}.kQ];
      case 'standard_rev_bi_uni',  this_v = [ r{it}.Vf;  r{it}.Vr;  r{it}.kA; r{it}.kB; r{it}.kP;];
      case 'standard_rev_bi_bi',   this_v = [ r{it}.Vf;  r{it}.Vr;  r{it}.kA; r{it}.kB; r{it}.kP; r{it}.kQ];
      case 'standard_rev_multiple', this_v = [ r{it}.Vf;  r{it}.Vr;  r{it}.kf; r{it}.kr;];
      case 'ready_made',            this_v = [ r{it}.Vf;  r{it}.Vr;  r{it}.KMf; r{it}.KMr; r{it}.KA; r{it}.KI];
    end

    if isfield(r{it},'Krep'), this_v = [this_v; r{it}.Krep]; end
    if isfield(r{it},'Kact'), this_v = [this_v; r{it}.Kact]; end

    rel_reactions(length(v)+1:length(v)+length(this_v)) = it;

    v = [v; this_v];

    for it2 = 1:length(kinetics.reactions{it}.parameters),
      names=[names; repmat({ ['R' num2str(it) ':' kinetics.reactions{it}.parameters{it2}]}, ...
                           kinetics.reactions{it}.sizes(it2),1) ];
    end
  end
end

for it =1:length(rel_reactions),
  relevant_reactions{it} = rel_reactions(it);
end
relevant_reactions =   relevant_reactions';