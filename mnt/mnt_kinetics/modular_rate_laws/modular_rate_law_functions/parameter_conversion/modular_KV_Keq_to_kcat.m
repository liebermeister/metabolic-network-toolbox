function [kcatplus, kcatminus] = modular_KV_Keq_to_kcat(N,kinetics,KV,Keq)

% [kcatplus, kcatminus] = modular_kcat_to_KV_Keq(N,kinetics,KV,Keq)

eval(default('KV','kinetics.KV','Keq','kinetics.Keq'));

all_KM              = ones(size(N'));
all_KM(find(N'~=0)) = kinetics.KM(find(N'~=0));
KMprod              = prod(all_KM.^(N'),2);

kcatplus  = KV .* sqrt(Keq.*KMprod);
kcatminus = KV ./ sqrt(Keq.*KMprod);
