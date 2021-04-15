function [kinetics,KV,Keq] = modular_kcat_to_KV_Keq(N,kinetics,kcatplus,kcatminus,Keq)

% [kinetics,KV,Keq] = modular_kcat_to_KV_Keq(N,kinetics,kcatplus,kcatminus,Keq)

eval(default('kcatplus','[]','kcatminus','[]','Keq','[]'));

if isempty(kcatminus),
  kcatminus = kinetics.Kcatr;
end

if isempty(kcatplus),
  kcatplus = kinetics.Kcatf;
end

all_KM              = ones(size(N'));
all_KM(find(N'~=0)) = kinetics.KM(find(N'~=0));
KMprod              = prod(all_KM.^(N'),2);
% note: KMprod corresponds to KM(products)/KM(substrates)

if length(kcatminus),
  Keq = kcatplus ./ kcatminus .* KMprod;
end

KV  = sqrt([kcatplus .* kcatminus]);

kinetics.KV  = KV;
kinetics.Keq = Keq;
