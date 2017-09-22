function [feasible,ind_non_orthogonal] = EBA_orth(sign_J,C)

% orthogonal = EBA_orth(sign_J,C)
%
% test flux sign vectors for thermodynamic feasibility
% function for comparing entire matrices -> see EBA_orth2
%
% sign_J sign vector (column)
% C matrix of cycle vectors (columns)
% feasible: flag stating that sign_J is sign-orthogonal to all cycles
% (i.e. all columns of C)

n_cycles       = size(C,2);
ind_common_support = find(column(abs(sign_J))' * abs(C));

CC              = C(:,ind_common_support);
A              = diag(sign_J)*CC;
orthogonal_list= (sum(A > 0) .* sum(A < 0))~=0;
feasible       = prod(double(orthogonal_list));

ind_non_orthogonal = [];
if ~feasible, ind_non_orthogonal = ind_common_support(find(orthogonal_list==0)); end
