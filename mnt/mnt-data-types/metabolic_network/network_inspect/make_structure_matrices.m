function [Mplus, Mminus, Wplus, Wminus, nm, nr, N_int,ind_M,ind_Wp,ind_Wm,ind_Mp,ind_Mm] = make_structure_matrices(N,W,ind_ext,h)

% [Mplus, Mminus, Wplus, Wminus, nm, nr, N_int,ind_M,ind_Wp,ind_Wm] = make_structure_matrices(N,W,ind_ext,h)

% nm       # metabolites
% nr       # reactions
%
% N        mn x nr Stoichiometric matric
% W        nr x nm Allosteric regulation matrix (+1 for activator, -1 for inhibitor) 
% ind_ext  vector of indices of external metabolites (numbers between 1 and nm)
% h        nr x 1 column vector of "cooperativity factors" (usually, they can all be set to 1)
%
% Mplus    nr x nm matrix of forward molecularities (ie |N|, but only for substrate elements)
% Mminus   nr x nm matrix of reverse molecularities (ie N, but only for product elements)
% Wplus    nr x nm matrix of activations (ie W, but only for activation elements)
% Wminus   nr x nm matrix of inhibitions (ie |W|, but only for inhibition elements)

eval(default('ind_first_order','[]','h','[]'));

N = sparse(N);
W = double(sparse(W));

[nm,nr] = size(N);

if isempty(h),
  h = ones(nr,1);
end

Mplus   = diag(h) * [abs(N').*(N'<0)];
Mminus  = diag(h) * [    N' .*(N'>0)];

Wplus   =     W .*(W>0);
Wminus  = abs(W).*(W<0);

N_int  = N(setdiff(1:nm,ind_ext),:);

ind_M  = find(Mplus+Mminus);
ind_Wp = find(Wplus);
ind_Wm = find(Wminus);
ind_Mp = find(Mplus);
ind_Mm = find(Mminus);
