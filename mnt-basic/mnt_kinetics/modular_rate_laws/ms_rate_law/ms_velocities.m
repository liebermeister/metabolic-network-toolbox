function  [v, v_plus, v_minus] = ms_velocities(N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h)

% [v, v_plus, v_minus] = ms_velocities(N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h)

[v, v_plus, v_minus] = modular_velocities('ms',N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h);