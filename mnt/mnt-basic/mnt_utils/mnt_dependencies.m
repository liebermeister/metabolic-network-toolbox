function mnt_dependencies(only_required)

% mnt_dependencies(only_required)
%
% Check whether the packages used by the Metabolic Network Toolbox are installed
% (see README.md for details)
%
% only_required: if true, check only matlab-utils and SBtab for Matlab,
%                and no optional packages (default: false)

if nargin < 1
  only_required = false;
end

if ~exist('matlab_utils_installed','file')
  warning('Please install the Matlab utility functions (https://github.com/liebermeister/matlab-utils) - They are required by the Metabolic Network Toolbox.');
end

if exist('package_version_check','file')
  package_version_check('The Matlab utility functions', 'matlab_utils_version', '1.0.0');
end

if ~only_required && ~exist('linprog','file')
  warning('Please install the Matlab Optimization Toolbox - Otherwise many flux analysis functions will not work.');
end

if ~only_required && ~exist('TranslateSBML','file')
  warning('Please install the Matlab interface of libSBML (https://sbml.org/software/libsbml/) - Otherwise the SBML import/export functions will not work.');
end

if ~exist('sbtab_version','file')
  warning('Please install SBtab for Matlab (https://github.com/liebermeister/sbtab-matlab) - Otherwise SBtab import and export will not work.');
end

if exist('package_version_check','file')
  package_version_check('SBtab for Matlab', 'sbtab_toolbox_version', '1.0.0');
end

if ~only_required && ~exist('sptensor','file')
  warning('Please install the Tensor Toolbox (https://www.tensortoolbox.org) - Otherwise second-order elasticities and response coefficients cannot be computed.');
end

if ~only_required && ~exist('CalculateFluxModes','file')
  warning('Please install efmtool (https://csb.ethz.ch/tools/software/efmtool.html) - Otherwise elementary flux modes cannot be computed.');
end
