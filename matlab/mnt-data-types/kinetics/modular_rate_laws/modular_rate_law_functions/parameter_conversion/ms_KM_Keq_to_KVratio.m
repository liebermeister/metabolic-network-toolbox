function kcat_ratio = ms_KM_Keq_to_KVratio(N,KM,Keq);

% kcat_ratio = ms_KM_Keq_to_KVratio(N,KM,Keq);

all_KM              = ones(size(N'));
all_KM(find(N'~=0)) = KM(find(N'~=0));
kcat_ratio          = Keq .* prod(all_KM.^(-N'),2);
