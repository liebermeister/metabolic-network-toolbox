Metabolic Network Toolbox for Matlab
====================================

The Metabolic Network Toolbox contains functions for constructing, simulating, and analysing
computational models of metabolism: stoichiometric network analysis, flux balance analysis,
thermodynamics, kinetic rate laws, simulation, metabolic control analysis, and graphics.

## Quick start

After the installation (see below), try this in Matlab:

```matlab
N = [-1 0; 1 -1; 0 1];                                % stoichiometric matrix
network = network_construct(N, [1; 1], [1 3]);        % all reversible; metabolites 1 and 3 external
network.kinetics = set_kinetics(network, 'cs');       % choose rate laws
[S, J] = network_steady_state(network, [2; 1; 1]);    % compute a steady state
```

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

For an overview of common tasks, see [GUIDE.md](GUIDE.md). The folder `mnt/demo` contains short examples
for single topics (see `help demo`).

## Requirements

- Matlab (tested with R2025b)
- [Matlab utility functions](https://github.com/liebermeister/matlab-utils) (required)
- [SBtab for Matlab](https://github.com/liebermeister/sbtab-matlab) (required for reading and writing
  SBtab files: kinetic data, state data, network layouts, and conversion between SBtab and networks)
- Optimization Toolbox (for linear and quadratic optimisation, used in many flux analysis functions).
  Optionally, Gurobi can be used instead of `linprog` and `quadprog` in the more demanding
  functions (see `help optim_solver` in matlab-utils).
- Statistics and Machine Learning Toolbox (only for `prob_of_order_relation`)

Some functions require further packages:

| Package | Needed for |
|---------|------------|
| [libSBML](https://sbml.org/software/libsbml/) (Matlab interface) | SBML import and export (`network_sbml_import`, `network_sbml_export`) |
| [Tensor Toolbox](https://www.tensortoolbox.org) | second-order elasticities and response coefficients |
| [efmtool](https://csb.ethz.ch/tools/software/efmtool.html) | elementary flux modes (`network_efmtool`, `eba_feasible`) |
| [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/) | conversion from and to COBRA models |
| [lasso.m](http://www.imm.dtu.dk/~hoe/files/lasso.m) | method `'lasso'` in `project_fluxes` (the function `lasso` of the Statistics Toolbox has a different form) |

## Contents

The folder `mnt` contains:

| Folder | Contents |
|---|---|
| `metabolic_network` | network structure: constructing, editing, inspecting |
| `kinetics` | rate laws: mass-action, convenience kinetics, modular rate laws, ... |
| `kinetic_data`, `state_data`, `omics_data` | kinetic data, metabolic states, omics data |
| `fluxes` | flux balance analysis, flux minimisation, flux sampling, energy balance analysis, max-min driving force |
| `thermodynamics` | driving forces and equilibrium constants |
| `simulation` | steady states and dynamic simulation |
| `control_analysis` | metabolic control analysis, spectral response coefficients (= frequency response coefficients; see below), models with uncertain parameters |
| `graphics` | network graphics |
| `io` | import and export: SBML, SBtab, COBRA (also models with a metabolic state, see below) |
| `model_and_state` | models with a reference state: loading, mapping data to a model |
| `annotation` | conversion between BiGG and KEGG identifiers (with a mapping table) |
| `data_orientation` | conversion of data tables to the reaction orientation of a model |
| `utils`, `demo` | constants, version numbers, demo script |

Networks, kinetics, and data are represented by Matlab structs; for details, type
`help metabolic_network`, `help kinetics`, `help kinetic_data`, and `help state_data`. Most folders contain a `Contents.m` file with an
overview (type, e.g., `help fluxes`).

Terminology: in function names and older help texts, "spectral response" and "spectral response
coefficients" mean frequency response and frequency response coefficients: the complex-valued
responses to sine-wave parameter perturbations, as functions of the frequency. The function
names keep the older term.

## Models with a metabolic state in SBML

SBML has no elements for flux data or equilibrium constants. To store a model together with a
metabolic state, MNT uses a customised use of standard SBML elements (SBML Level 3 Version 2,
functions `sbml_network_state_save` and `sbml_network_state_load`; the files are valid SBML and
can be read by any SBML tool):

- concentrations: initial concentrations of the species (default unit mM);
- fluxes, equilibrium constants and reaction Gibbs free energies: local parameters `reference_flux`
  (SBO:0000612, default mM/s), `equilibrium_constant` (SBO:0000281, mM to the power of the sum of
  the stoichiometric coefficients) and `reaction_gibbs_free_energy` (SBO:0000617, kJ/mol) in the
  kinetic law of each reaction, with unit definitions and notes explaining their meaning;
- for networks without kinetics, the kinetic law is the constant rate
  `reference_flux * flux_reference_volume` (1 litre), so a simulation stays in the reference state.

Compartments, database identifiers (KEGG, BiGG, ChEBI, as identifiers.org annotations) and the
original metabolite and reaction ids are kept. The same states are stored in SBtab files by
`sbtab_network_state_save` (tables `MetaboliteConcentration`, `Flux`, `EquilibriumConstant`,
`ReactionGibbsFreeEnergy`); `sbtab_to_sbml` and `sbml_to_sbtab` convert between the two formats,
and `mnt_load_model_and_state` loads a model with its state from a .mat, SBtab or SBML file.
Details: `help sbml`.

## Corrections (November 2026)

A check of the modular rate law formulae (Liebermeister et al., Bioinformatics 2010) against the code
led to these corrections (details in `CHANGELOG.md` and in the help texts of the functions):

- `compute_modular_elasticities`: the rate laws `'ds'` (DM) and `'fd'` (FM) stopped with an error (undefined
  variables `theta_plus`, `theta_minus`, and a wrong logarithm, log(1/(alpha-1)) instead of log(1/alpha-1)).
  The other rate laws (`'cs'`, `'ms'`, `'rp'`) were correct and are unchanged.
- `compute_modular_response_second`: second-order response coefficients of a target function were correct only if
  the target is a single concentration or flux (the usual case). For linear combinations of concentrations and fluxes,
  the term -R^z R^z was replaced by a weighted sum of -R^y R^y. The new argument `target_type` ('linear', default, or
  'loglinear') states how the target is defined.
- The first- and second-order elasticities of all five rate laws, and the second-order response coefficients, are now
  tested against numerical derivatives (`tests/test_modular_elasticities.m`).

The functions for control and response coefficients (`control_coefficients`, `response_coefficients`,
`response_coefficients_sparse`, `norm_response_coefficients`) and the conversion between scaled and unscaled
second-order elasticities were checked and found to be correct.

## Tests

From the repository folder, run

```matlab
results = runtests('tests')
```

## Versions

To see the installed version, type `mnt_version`. Older versions can be downloaded from the [releases page](https://github.com/liebermeister/metabolic-network-toolbox/releases) of this repository, and changes between versions are listed in `CHANGELOG.md`.

Other packages depend on this toolbox, so functions are not simply renamed or removed: an old function
name remains available for at least one MINOR version as a thin wrapper that calls the new function and
gives a warning. The name is removed only in the next MAJOR version. Functions marked as unfinished in
their help text are not covered by this rule.

## Contributing

Conventions for code and help texts are described in `CONTRIBUTING.md`.

## License
This package is released under the [GNU General Public License](LICENSE), version 2 or later.

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.

## How to cite
Citation information is provided in the file `CITATION.cff`.
