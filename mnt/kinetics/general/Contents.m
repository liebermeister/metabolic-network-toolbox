% General functions for rate laws
%
%   kinetics_print                    - Print a kinetics structure
%   kinetics_string                   - Represent kinetics of the nth reaction as a string
%   kinetics_to_strings               - Convert a kinetics structure to kinetics of type 'kinetic_strings'
%   network_execute_initialAssignment - Execute network.kinetics.initialAssignment and update parameter values and compartment sizes
%   network_get_kinetics_strings      - Get kinetics strings and names of the rate laws for a network with modular rate laws
%   network_get_parameter_values
%   network_print_kinetics            - Display the rate laws of a network
%   network_set_enzyme_levels         - Set the enzyme levels in network.kinetics
%   network_velocities                - Compute reaction rates for given metabolite concentrations
%   parameters2vector                 - Convert kinetic parameters to a vector
%   set_kinetics                      - Construct the field 'kinetics' of a network
%   set_standard_kinetics             - Construct a kinetics field for a metabolic network
%   standard_kinetics2vector          - Extract parameter vector from kinetics of type 'standard'
%   standard_print                    - Display kinetic rate laws
%   vector2parameters                 - Convert a parameter vector to a kinetics structure
%   vector2standard_kinetics          - Construct kinetics structure from a parameter vector (for kinetics of type 'standard')
%   vectorchange2standard_kinetics    - Construct kinetics structure from a parameter vector (for kinetics of type 'standard')
