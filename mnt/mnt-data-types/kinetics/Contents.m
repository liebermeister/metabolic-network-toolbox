% KINETICS Data structure 'kinetics' for rate laws in metabolic networks
%
% The rate laws of a network are stored in the field network.kinetics (see 'help metabolic_network').
% To construct this field, use set_kinetics. To compute reaction rates, use network_velocities.
% To display the rate laws, use network_print_kinetics.
%
% The field kinetics.type shows the type of rate law. With m metabolites and n reactions:
%
% Modular rate laws: type 'cs', 'ms', 'ds', 'rp', or 'fd'
%   (common saturable, simultaneous binding, direct binding, power law, and force-dependent)
%   u      (n x 1 vector)    enzyme levels
%   c      (m x 1 vector)    metabolite levels
%   KM     (n x m matrix)    Michaelis constants
%   KA     (n x m matrix)    activation constants
%   KI     (n x m matrix)    inhibition constants
%   KV     (n x 1 vector)    velocity constants
%   Keq    (n x 1 vector)    equilibrium constants
%   Kcatf  (n x 1 vector)    forward turnover rates
%   Kcatr  (n x 1 vector)    reverse turnover rates
%   h      (n x 1 vector)    reaction cooperativities
%
% Type 'convenience': convenience kinetics
%   r      (n x 1 vector)    geometric mean of forward and backward turnover rates
%   g      (m x 1 vector)    energy constants, g = exp((G_formation-G_mean)/RT * scaling factor)
%                            (for the scaling factor, see convert_G_to_kG and convert_kG_to_G)
%   KM, KA, KI (n x m matrices)   as above
%   E      (n x 1 vector)    enzyme levels
%   S      (m x 1 vector)    metabolite levels
%
% Type 'mass-action': mass-action kinetics
%   k_fwd      (n x 1 vector)  forward rate constants
%   k_bwd      (n x 1 vector)  backward rate constants
%   exponents  (optional)      for reactions of type 2A -> B
%
% Type 'numeric': rate laws given by a Matlab function
%   velocity_function          function handle: v = my_function(s,p) or v = my_function(s,p,t),
%                              with column vectors s and v
%   parameters                 struct with the parameters p
%
% Type 'kinetic_strings': rate laws given as strings
%   reactions                  (n x 1 cell) one struct for each reaction
%   parameters                 (cell) global parameter names
%   parameter_values           (vector) global parameter values
%
% Type 'standard': each reaction has its own type of rate law
%   reactions                  (n x 1 cell) one struct for each reaction; its field 'type'
%                              gives the rate law, for instance 'mass-action' or
%                              'michaelis_menten_rev'
%
% Conversion between parameters and vectors: parameters2vector, vector2parameters
