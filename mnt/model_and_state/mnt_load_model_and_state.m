function [network, v, c, keq, theta, networkG] = mnt_load_model_and_state(model_state_file, filetype)

% [network, v, c, keq, theta, networkG] = mnt_load_model_and_state(model_state_file, filetype)
%
% Load a model (network struct) and its reference state from a .mat, SBtab or SBML file
%
% model_state_file  file name; the file type is recognised by the extension:
%                     .mat          Matlab file with variables network, v, c, keq (and optionally
%                                   networkG; a field of a struct is not enough)
%                     .tsv, .csv    SBtab file (sbtab_network_state_load)
%                     .xml, .sbml   SBML file  (sbml_network_state_load)
% filetype          (optional) 'matlab', 'sbtab' or 'sbml', for files with other extensions
%
% Outputs
%   network   network struct
%   v         fluxes
%   c         metabolite concentrations
%   keq       equilibrium constants; for SBtab and SBML files without equilibrium constants,
%             computed from the reaction Gibbs free energies
%   theta     thermodynamic forces -dG/RT: theta = log(keq) - network.N' * log(c)
%             (fluxes v must be conformal with theta)
%   networkG  network struct with graphics information (only in .mat files; otherwise [])

eval(default('filetype','[]'));

if isempty(filetype)
  [~, ~, ext] = fileparts(model_state_file);
  switch lower(ext)
    case '.mat',            filetype = 'matlab';
    case {'.tsv','.csv'},   filetype = 'sbtab';
    case {'.xml','.sbml'},  filetype = 'sbml';
    otherwise
      error('Unknown file extension "%s"; please give the file type (''matlab'', ''sbtab'' or ''sbml'')', ext);
  end
end

networkG = [];
dmu      = [];

switch filetype
  case 'matlab'
    X        = load(model_state_file);
    network  = X.network;
    v        = X.v;
    c        = X.c;
    keq      = X.keq;
    if isfield(X,'networkG'), networkG = X.networkG; end
  case 'sbtab'
    [network, c, v, dmu, keq] = sbtab_network_state_load(model_state_file);
  case 'sbml'
    [network, c, v, dmu, keq] = sbml_network_state_load(model_state_file);
  otherwise
    error('Unknown file type "%s" (allowed: ''matlab'', ''sbtab'', ''sbml'')', filetype);
end

if isempty(keq)
  if isempty(dmu)
    error('The file contains neither equilibrium constants nor reaction Gibbs free energies');
  end
  dmu0 = dmu - RT * network.N' * log(c);
  keq  = exp(-dmu0/RT);
end

theta = log(keq) - network.N' * log(c);
