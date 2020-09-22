function [kcatplus,kcatminus] = ms_compute_Kcat(N,KM,KV,Keq,h);

% [kcatplus,kcatminus] = ms_compute_Kcat(N,KM,KV,Keq);
%
% Let
% nm       # metabolites
% nr       # reactions
%  
% Input:  
% N        mn x nr Stoichiometric matric
% KM       nr x nm sparse matrix of KM values (Michaelis-Menten)
% KV       nr x 1 column vector of KV values (velocity constants, defined as sqrt(Kcat_forward * Kcat_reverse)
% Keq      nr x 1 column vector of equilibrium constants
% h        nr x 1 column vector of "cooperativity factors" (usually, they can all be set to 1)
%
% Output: 
% kcatplus  nr x 1 column vector of forward kcat values
% kcatminus nr x 1 column vector of reverse kcat values
  
%  
if ~exist('h','var'), h = ones(size(KV)); end
  
[log_kcatplus,log_kcatminus] = ms_compute_log_Kcat(N,KM,KV,Keq,h);

kcatplus  = exp(log_kcatplus);
kcatminus = exp(log_kcatminus);

