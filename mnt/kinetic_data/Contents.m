% KINETIC_DATA Data structure 'kinetic_data' for kinetic and thermodynamic data of a network
%
% The struct contains one field for each quantity, named by its symbol (for instance Keq, KM, or
% dmu0). To construct it, use kinetic_data_load, or kinetic_data_construct for an empty structure.
%
% Each quantity field contains arrays of the same size:
%   scaling            'Additive' or 'Multiplicative': type of scale for this quantity
%   median, mean, std  values on the original scale
%   lower, upper       lower and upper bounds
%   mean_ln, std_ln    mean and standard deviation of the logarithm
%   geomstd            geometric standard deviation
%   lower_ln, upper_ln bounds of the logarithm
%
% Functions:
%   kinetic_data_load         read data from an SBtab file
%   kinetic_data_save         save data to an SBtab file
%   kinetic_data_print        display data
%   kinetic_data_complete     insert bounds and pseudo values
%   kinetic_data_subselect    select data for a submodel
%   kinetic_data_to_kinetics  convert data into a 'kinetics' structure
%   kinetics_to_kinetic_data  convert a 'kinetics' structure into data
%   mnt_prior_file            file name of the default prior table
%   parameter_balancing_prior read a prior table
