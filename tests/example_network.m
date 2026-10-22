function network = example_network()

% network = example_network()
%
% Network from demo_mnt, used in the tests:
%
%  S1             S5
%     \         /
%       S3 - S4
%     /         \
%  S2             S6
%
% S1, S2, S5, and S6 are external metabolites; all reactions are reversible

N = [ -1  0  0  0  0; ...
       0 -1  0  0  0; ...
       1  1 -1  0  0; ...
       0  0  1 -1 -1; ...
       0  0  0  1  0; ...
       0  0  0  0  1];

network = network_construct(N, [1 1 1 1 1]', [1 2 5 6]');
