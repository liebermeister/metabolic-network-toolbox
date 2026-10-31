function keq = ms_KM_KVratio_to_Keq(N,KM,kcat_ratio);

% keq = ms_KM_KVratio_to_Keq(N,KM,kcat_ratio);

all_KM              = ones(size(N'));
all_KM(find(N'~=0)) = KM(find(N'~=0));
keq                 = kcat_ratio .* prod(all_KM.^(N'),2);
