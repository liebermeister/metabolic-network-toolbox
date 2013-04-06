%MNT_DEMO - Demo for Metabolic Network Toolbox

data_dir = [fileparts(which(mfilename)) '/data'];

echo on;
clc
% We consider the following simple network
%
% S1             S5
%    \         /
%      S3 - S4
%    /         \
% S2             S6
%
% (where Sn denotes metabolites) to illustrate some
% of the functions of the Metabolic Network Toolbox.
%
% Please press any key to continue

pause
clc
% We define the network by its stoichiometric matrix N,
% the reversibility of reactions (bit string) 
% and the external metabolites (indices)
%

N            = [ -1 0 0 0 0; ...
		 0 -1 0 0 0; ...
		 1 1 -1 0 0; ...
		 0 0 1 -1 -1; ...
		 0 0 0 1 0; ...
		 0 0 0 0 1];

reversible   = [1 1 1 1 1]';

external_ind = [1 2 5 6 ]';

%
% Press any key to continue
pause

% All information about the network is then stored in a matlab structure

network      = network_construct(N,reversible,external_ind)
% The fields 'kinetics' and 'graphics_par' contain information
% about reaction kinetics and the graphical layout (set here by default).

% Press any key to continue
pause
clc

% Now we analyse the network structure
% Information about stationary fluxes and conservation
% relations can be obtained from the kerbel matrix K
% and the link matrix L (for the stoichimetric matrix 
% of internal metabolites)
[ K, L ] = analyse_N(N(find(network.external==0),:))
% Press any key to continue
pause
clc
% Until here, we only considered the network structure.
% Now we define a kinetics.
% There are different ways to specify a kinetics (see 'kinetics_structure'). 
% Here we choose a common modular rate law kinetics for all reactions.
%
network.kinetics = set_kinetics(network,'cs');
% This updates the field 'kinetics', which contains all kinetic data.
%
% Press any key to continue
pause
clc
% Now we can run dynamical simulations.
% First, we compute a stationary state, starting from an initial guess
% of the concentrations ...

s = [1 1 0 0 0 0]';

[S, J] = network_steady_state(network,s);

% ... and plot it on the network graph. 
%

network = netgraph_make_graph(network);
netgraph_draw(network,{'metvalues',S,'arrowvalues',J,'arrowstyle','fluxes'});

% The size (and color) of octagons 
% and squares represents the values of concentrations and fluxes.
%
% Press any key to continue
pause
clc
% Likewise, we can compute properties from metabolic control analysis,
% for instance, elasticity coefficients and control coefficients
% (for more information: type 'help control_theory')
%

epsilon       = elasticities(network,S);
[ C_J, C_S ]  = control_coefficients(network.N,epsilon,network.external);

% Here we plot the control coefficients on the 4th metabolite 

netgraph_draw(network,'actvalues',C_S(4,:)');

% Press any key to continue
pause
clc
% Metabolic control theory has proven theorems that relate the 
% control coefficients to the network structure (described by 
% K and L) and the reaction elasticities.
%
% The theorems for non-normalised control coefficients read
%
%Summation theorems
% CJ K         = K      
% CS K         = 0
%
%Connectivity theorems
% CJ epsilon L = 0
% CS epsilon L = - L
%
%Let us check them:
%
% CJ K = K

C_J*K
K
pause
% CS K = 0
C_S(find(network.external==0),:)*K
pause
% CJ epsilon L = 0
C_J(find(network.external==0),:)*epsilon(:,find(network.external==0))*L
pause
% CS epsilon L = - L
C_S(find(network.external==0),:)*epsilon(:,find(network.external==0))*L
L
pause
clc
% We compute a time course with initial conditions given by s:

[t,s_t,s_int] = network_integrate(network,s,5);
close
plot(t,s_t)

% Press any key to continue
pause

% We can also show the time course as a movie
M = netgraph_movie(network,t,s_t);
movie(M,1);
% Press any key to continue
pause
close

% Enjoy this toolbox!
echo off;
