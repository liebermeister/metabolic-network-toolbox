function [Kplus,Kminus] = ms_compute_Kcat(N,KM,KV,Keq);

% [Kplus,Kminus] = ms_compute_Kcat(N,KM,KV,Keq);

[log_Kplus,log_Kminus] = ms_compute_log_Kcat(N,KM,KV,Keq);
Kplus                  = exp(log_Kplus);
Kminus                 = exp(log_Kminus);

