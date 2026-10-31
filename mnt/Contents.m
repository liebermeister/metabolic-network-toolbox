% Metabolic Network Toolbox
% Version 1.3.0 31-Oct-2026
%
% Functions for constructing, simulating, and analysing computational models of metabolism.
% Models are represented by Matlab structs. For a first example, run demo_mnt.
% For a list of common tasks, see GUIDE.md in the main folder.
%
% Data structures
%   metabolic_network   - network structure: constructing, editing, inspecting ('help metabolic_network')
%   kinetics            - rate laws ('help kinetics')
%   kinetic_data        - kinetic and thermodynamic data ('help kinetic_data')
%   state_data          - metabolite levels, fluxes, and enzyme levels ('help state_data')
%   omics_data          - omics data ('help biotable')
%
% Analysis
%   fluxes              - flux balance analysis, thermodynamic constraints, flux sampling
%   thermodynamics      - driving forces and equilibrium constants
%   simulation          - steady states and dynamic simulation
%   control_analysis    - metabolic control analysis, spectral response, uncertain parameters
%
% Other
%   graphics            - network graphics ('help netgraph_draw')
%   io                  - import and export: SBML, SBtab, COBRA
%   utils               - constants and version numbers
%   demo                - demo script
%
% Other packages required: see README.md
%
% Copyright (C) 2006-2026
% Wolfram Liebermeister <wolfram.liebermeister@gmail.com>
