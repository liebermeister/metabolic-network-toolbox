function [r,q] = network_issplit(network)

% [split_flag,split_matrix] = network_issplit(network)
%
% split_matrix : matrix relating irreversible semireactions to reversible  reactions

r = 0;
q = eye(length(network.actions));

if isfield(network,'split_reactions'),
  r = network.split_reactions;
  q = network.split_matrix;
end