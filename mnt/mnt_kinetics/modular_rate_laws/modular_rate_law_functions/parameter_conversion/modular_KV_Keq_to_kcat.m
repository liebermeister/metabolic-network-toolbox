function [kcatplus, kcatminus] = modular_KV_Keq_to_kcat(N,kinetics,KV,Keq,KM,h)

% [kcatplus, kcatminus] = modular_kcat_to_KV_Keq(N,kinetics,KV,Keq,KM)

eval(default('KV','kinetics.KV','Keq','kinetics.Keq','KM','kinetics.KM','h','ones(size(KV))'));

[kcatplus,kcatminus] = ms_compute_Kcat(N,KM,KV,Keq);

% DEPRECATED
% all_KM              = ones(size(N'));
% all_KM(find(N'~=0)) = kinetics.KM(find(N'~=0));
% prod_KM             = prod(all_KM.^(N'),2);
% kcatplus            = KV .* [sqrt(Keq ./ prod_KM) .^ h];
% kcatminus           = KV ./ [sqrt(Keq ./ prod_KM) .^ h];
