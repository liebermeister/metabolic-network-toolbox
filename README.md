Metabolic Network Toolbox for Matlab
====================================

The Metabolic Network Toolbox contains functions for constructing, simulating, and analysing
computational models of metabolism: stoichiometric network analysis, flux balance analysis,
thermodynamics, kinetic rate laws, simulation, metabolic control analysis, and graphics.

## Installation

Download or clone the repository and add the `mnt` folder, with its subfolders, to your Matlab path:

```matlab
addpath(genpath('/path/to/metabolic-network-toolbox/mnt'));
savepath;   % optional: keep the folders on the path in future sessions
```

Then install the required packages (see below) and run the demo:

```matlab
demo_mnt
```

To check which of the other packages are installed, run `mnt_dependencies`.

## Requirements

- Matlab (tested with R2025b)
- [Matlab utility functions](https://github.com/liebermeister/matlab-utils) (required)
- Optimization Toolbox (for linear and quadratic optimisation, used in many flux analysis functions).
  Optionally, Gurobi can be used instead of `linprog` and `quadprog` in the more demanding
  functions (see `help optim_solver` in matlab-utils).

Some functions require further packages:

| Package | Needed for |
|---------|------------|
| [libSBML](https://sbml.org/software/libsbml/) (Matlab interface) | SBML import and export (`network_sbml_import`, `network_sbml_export`) |
| [SBtab for Matlab](https://github.com/liebermeister/sbtab-matlab) | reading and writing SBtab files (kinetic data, state data, network layouts) |
| [Tensor Toolbox](https://www.tensortoolbox.org) | second-order elasticities and response coefficients |
| [efmtool](https://csb.ethz.ch/tools/software/efmtool.html) | elementary flux modes (`network_efmtool`, `eba_feasible`) |
| [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/) | conversion from and to COBRA models |

## Contents

The functions are organised in three folders:

- **`mnt/mnt-data-types`**: data structures for metabolic networks, kinetics (mass-action,
  convenience kinetics, modular rate laws, ...), kinetic data, and metabolic state data;
  import and export (SBML, SBtab, COBRA)
- **`mnt/mnt-basic`**: flux analysis (FBA, flux minimisation, flux sampling, energy balance
  analysis, max-min driving force), thermodynamics, simulation and steady states, graphics,
  conversion between SBtab and networks, omics data, SBML and COBRA (`mnt_sbtab_conversion`)
- **`mnt/mnt-control-analysis`**: metabolic control analysis, spectral response coefficients,
  models with uncertain parameters

Networks and their kinetics are represented by Matlab structs; for details, type
`help metabolic_network` and `help kinetics`. Each folder contains a `Contents.m` file with an
overview (type, e.g., `help mnt-basic`).
HTML documentation is available in the folder `doc`.

## Tests

From the repository folder, run

```matlab
results = runtests('tests')
```

## Versions

To see the installed version, type `mnt_version`. Older versions can be downloaded from the [releases page](https://github.com/liebermeister/metabolic-network-toolbox/releases) of this repository, and changes between versions are listed in `CHANGELOG.md`.

## License
This package is released under the [GNU General Public License](LICENSE), version 2 or later.

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.
