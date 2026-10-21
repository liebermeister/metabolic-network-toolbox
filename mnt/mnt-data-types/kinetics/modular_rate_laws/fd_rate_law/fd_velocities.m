function  [v, v_plus, v_minus] = fd_velocities(N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h)

[v, v_plus, v_minus] = modular_velocities('fd',N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,h);