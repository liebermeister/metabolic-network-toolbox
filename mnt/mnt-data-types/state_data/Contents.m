% STATE_DATA Data structure 'state_data' for metabolic states
%
% The struct contains metabolite levels, fluxes, and enzyme levels for one or several
% metabolic states, mapped onto a network (see 'help metabolic_network').
% In all matrices, the columns refer to the states:
%
%   samples                     names of the states
%   metabolite_data.Mean, .Std  (m x #states matrices) metabolite levels
%   flux_data.Mean, .Std        (n x #states matrices) fluxes
%   enzyme_data.Mean, .Std      (n x #states matrices) enzyme levels
%
% Functions:
%   load_state_data_for_network  load metabolite, flux, or enzyme data from an SBtab file
%   state_data_save              save state data to an SBtab file
%   flux_data_load               load one flux distribution
%   data_to_state_data           convert data from model balancing into this structure
%   state_data_load              (outdated) load data for one state from separate files
