function mnt_dependencies()

% mnt_dependencies()
%
% Check whether the packages used by the Metabolic Network Toolbox are installed
% (see README.md for details)

if ~exist('matlab_utils_installed','file')
  warning('Please install the Matlab utility functions (https://github.com/liebermeister/matlab-utils) - They are required by the Metabolic Network Toolbox.');
end

if exist('package_version_check','file')
  package_version_check('The Matlab utility functions', 'matlab_utils_version', '1.0.0');
end

if ~exist('linprog','file')
  warning('Please install the Matlab Optimization Toolbox - Otherwise many flux analysis functions will not work.');
end

if ~exist('TranslateSBML','file')
  warning('Please install the Matlab interface of libSBML (https://sbml.org/software/libsbml/) - Otherwise the SBML import/export functions will not work.');
end

if ~exist('sbtab_version','file')
  warning('Please install SBtab for Matlab (https://github.com/liebermeister/sbtab-matlab) - Otherwise SBtab import and export will not work.');
end

if exist('package_version_check','file')
  package_version_check('SBtab for Matlab', 'sbtab_toolbox_version', '1.0.0');
end

if ~exist('sptensor','file')
  warning('Please install the Tensor Toolbox (https://www.tensortoolbox.org) - Otherwise second-order elasticities and response coefficients cannot be computed.');
end

if ~exist('CalculateFluxModes','file')
  warning('Please install efmtool (https://csb.ethz.ch/tools/software/efmtool.html) - Otherwise elementary flux modes cannot be computed.');
end
