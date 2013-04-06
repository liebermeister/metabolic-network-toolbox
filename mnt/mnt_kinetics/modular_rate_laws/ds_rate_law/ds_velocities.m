function  [v, v_plus, v_minus] = ds_velocities(N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h)

% [v, v_plus, v_minus] = ds_velocities(N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h)

[v, v_plus, v_minus] = modular_velocities('ds',N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h);