# How do I ...?

A short guide to the most common tasks. Type `help` and a function name for details. The short
examples in the folder `mnt/demo` can be opened with `edit`, for instance `edit mnt_example_fluxes`.

The main data structures are described in `help metabolic_network`, `help kinetics`,
`help kinetic_data`, and `help state_data`.

## First steps

```matlab
N = [-1 0; 1 -1; 0 1];                                % stoichiometric matrix
network = network_construct(N, [1; 1], [1 3]);        % all reversible; metabolites 1 and 3 external
network.kinetics = set_kinetics(network, 'cs');       % choose rate laws
[S, J] = network_steady_state(network, [2; 1; 1]);    % compute a steady state
```

## Build and change a network

Example: `mnt_example_network`

| Task | Function |
|---|---|
| Build a network from a stoichiometric matrix | `network_construct` |
| Build a network from a list of reaction formulae | `network_build_from_sum_formulae_list` |
| Build a network from formulae in an SBtab file | `network_build_from_sum_formulae` |
| Set metabolites to be external | `network_set_external` |
| Choose a part of a network | `network_choose`, `network_subnetwork` |
| Remove metabolites or reactions | `network_remove`, `network_remove_action` |
| Merge metabolites | `network_merge_metabolites` |

## Inspect a network

| Task | Function |
|---|---|
| Show the reaction formulae | `network_print_formulae` |
| Find metabolites or reactions by name | `network_find_metabolites`, `network_find_reactions` |
| Count reactions and metabolites | `network_numbers` |
| Compute kernel matrix and link matrix | `network_analyse`, `analyse_N` |
| Check a network structure | `network_check`, `network_check_consistency` |

## Rate laws

Example: `mnt_example_kinetics`

| Task | Function |
|---|---|
| Choose rate laws for a network | `set_kinetics` |
| Show the rate laws | `network_print_kinetics` |
| Compute reaction rates | `network_velocities` |
| Convert parameters to a vector, and back | `parameters2vector`, `vector2parameters` |

## Simulation

Example: `mnt_example_kinetics`

| Task | Function |
|---|---|
| Compute a steady state | `network_steady_state` |
| Simulate a time course | `network_integrate` |
| Simulate with an oscillating parameter | `network_integrate_kin` |

## Fluxes and thermodynamics

Example: `mnt_example_fluxes`

| Task | Function |
|---|---|
| Flux balance analysis | `fba`, with options from `fba_default_options` |
| Minimise the sum of absolute fluxes | `pmf` |
| Project fluxes to the stationary subspace | `project_fluxes` |
| Sample fluxes or chemical potentials | `sample_feasible_v`, `sample_feasible_mu` |
| Test whether fluxes are thermodynamically feasible | `eba_feasible_lp`, `eba_feasible` |
| Compute elementary flux modes (requires efmtool) | `network_efmtool` |
| Max-min driving force | `mdf` |

Linear and quadratic problems use `linprog` and `quadprog`. To use Gurobi instead, see
`help optim_solver` in matlab-utils.

## Metabolic control analysis

Example: `mnt_example_control_analysis`

| Task | Function |
|---|---|
| Compute the main quantities with one command | `basic_control_analysis` |
| Compute elasticities | `elasticities` |
| Compute control coefficients | `control_coefficients` |
| Compute response coefficients | `response_coefficients` |
| Compute spectral response coefficients | `spectral_response` |

## Graphics

Example: `mnt_example_graphics`

| Task | Function |
|---|---|
| Define a layout | `netgraph_make_graph` |
| Draw a network | `netgraph_draw` (all display settings are listed in its help text) |
| Show concentrations and fluxes | `netgraph_concentrations` |
| Move elements with the mouse | `netgraph_move` |
| Read or save a layout | `netgraph_read_positions`, `netgraph_print_positions` |
| Make a movie | `netgraph_movie`, `netgraph_flux_movie` |

## Import, export, and data

Example: `mnt_example_import_export`

| Task | Function |
|---|---|
| Load or save a model in SBtab format | `sbtab_to_network`, `network_to_sbtab` |
| Load a model with kinetic data and state data | `sbtab_load_network_model_and_data` |
| Save a model with concentrations and fluxes | `sbtab_network_state_save` |
| Load or save a model in SBML format (requires libSBML) | `network_sbml_import`, `network_sbml_export` |
| Convert from and to COBRA models | `cobra_to_network`, `network_to_cobra` |
| Load kinetic data | `kinetic_data_load` |
| Load metabolite, flux, or enzyme data | `load_state_data_for_network` |
| Work with omics data | see `help biotable` |
