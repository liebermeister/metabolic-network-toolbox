# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## Not released yet

- `mnt_BASEDIR` renamed to `mnt_basedir`; it now returns the main directory of the repository, without a final file separator. The old name cannot remain as a wrapper, because the file names differ only by case (it had no callers)
- New function folders have a `Contents.m` file; `README.md` has a Quick start section
- `sample_feasible_mu`: default method 'sample' now samples uniformly by hit-and-run (was: 'extreme_points'); methods 'sample' and 'centre' returned wrong vectors; works with rank-deficient stoichiometric matrices; method 'centre' can also be called 'central'
- `stm_default_options`: new field `mu_sampling_method`; removed unused fields `n_flux_samples`, `n_saturation_samples`, `method_flux_sampling`
- `modular_velocities`, `modular_KV_Keq_to_kcat` and `modular_fast_parameters` pass the cooperativity `h` to `ms_compute_Kcat`
- `sample_feasible_v`: predefined non-zero fluxes and bounds `v_min > 0` or `v_max < 0` are respected
- `compute_modular_elasticities`: reactions at equilibrium (v = 0, A = 0) get their elasticities from the one-way rates
- `compute_modular_response_second`: corrected second-order response coefficients (scaling, factor 1/z), `zv_scaled` default, and external metabolites
- `control_coefficients`: error for unstable states, or nan values with the new argument `on_unstable = 'nan'`; new output `stable`; eigenvalues that are positive by numerical error only are corrected, with a tolerance relative to the size of M
- `response_coefficients_sparse`: fixed crash in small networks (dense intermediate results)
- `stm_default_options`: `dmu_limit` lowered from 100 to 50 kJ/mol; new fields `flag_control` and `n_max_unstable`
- New folder structure, with fewer levels and one folder for import and export (see `README.md`). Function names do not change, so calls are not affected, as long as the toolbox is added to the path with its subfolders
- `help fba`, `help eba`, and `help mdf` now show the functions: their folders were renamed
- `GUIDE.md` with an overview of common tasks; short examples for single topics in the folder `demo`, which are run by the tests
- `CONTRIBUTING.md` with conventions for code and help texts
- Reference for the data structures (`help metabolic_network`, `help kinetics`, `help kinetic_data`, `help state_data`)
- Help texts of the main functions now have a one-line description and links to related functions, some with examples
- Network graphics (folder `graphics`):
  - `netgraph_draw`: **figures may look different**: with `showsign = 1` (default), `metvaluesmax`/`metvaluesmin` (and the same for reactions, edges, and regulation) now set the scale of symbols and colours; before, values were always scaled by their own maximum, so frames of `netgraph_movie` did not share a scale. The default colorbar now shows the symmetric range of the colours, [-M, M]
  - `netgraph_draw`: curved lines (`straightlines = 0`) work again; edges are drawn without arrows too, and arrows follow the curves
  - `netgraph_draw`: `single_arrow = 1` also shows arrows for negative fluxes
  - `netgraph_move` no longer stores `arrowstyle = 'none'` in the network (before, arrows disappeared from later drawings)
  - `netgraph_print_positions`: policies 'replace elements' (default) and 'add nonexistent' work if model and file contain different elements; if the file does not exist yet, it is created
  - `netgraph_find_x`: fixed positions work in networks with more than 255 nodes and stay exactly in place; the state of the random number generator is restored
  - `netgraph_split_metabolites` works for networks without `regulation_matrix`
  - Movie functions: time bar in `netgraph_flux_movie` fixed; background colours are used also if they contain zeros, and are interpolated in `netgraph_flux_movie`; `netgraph_flux_arrows_movie` has default options; EPS files show the current figure (before: figure 1)
  - `netgraph_read_positions` without the SBtab package works again
  - `help graphics` lists all functions; test file `tests/test_graphics.m`
- Network graphics: layout (**figures look different**):
  - `netgraph_draw` and `netgraph_concentrations` work in subplots: axes created by `subplot` (or placed by the user) keep their position; before, they were moved to fill the figure, and `netgraph_concentrations` deleted the other subplots. A new axes object is still enlarged to fill the figure
  - Axis limits contain all nodes, symbols, and labels (before: `axis tight`, which cut off labels); the option `figure_axis` now sets the limits
  - Name labels are placed to the right of the symbols, vertically centred (option `text_offset`, default now depends on `squaresize`)
  - Sizes (`squaresize`, `arrowsize`, text offsets, `circle_shift`, `shade_long_lines`) refer to a layout of size 1 and are scaled to the actual layout size, so layouts in other units (e.g. pixels) look the same; new option `scale_sizes` (default 1)
  - Colorbar outside the network, below the axes it belongs to, with automatic ticks; colormaps are set per axes; default `colorbar_fontsize` 10 (was 14)
  - Default colours: dark grey lines and node dots (were blue), grey direction arrows in `netgraph_draw`, red flux arrows in `netgraph_concentrations`; white figure background
  - `netgraph_make_graph`: default `arrowsize` is smaller in larger networks; metabolite names keep 'S-' and letters before hyphens (e.g. 'ACETYL-COA' was shown as 'ACETYCOA'); the stereo prefixes D-, L-, (R)-, (S)- are still omitted
  - `netgraph_read_positions`: elements missing in the layout file are placed near their neighbours in the network (before: all at the same point)
  - `netgraph_draw`: unknown fluxes (nan values in `arrowvalues`) are shown by grey diamonds (before: nothing was drawn)
  - `netgraph_make_graph`: the argument `table_positions` (layout file or SBtab object) is used; before, it was ignored
- Network graphics: options
  - All drawing functions take graphics options in the same way: as a struct (called `gp` in all help texts; before also `options`, `goptions`, `pars`) or as name-value pairs at the end of the argument list, e.g. `netgraph_concentrations(network, c, v, 1, 'colorbar', 1)`. In `netgraph_concentrations`, `netgraph_flux_movie`, and `netgraph_flux_arrows_movie`, options can also be given in place of `text_flag`. Existing calls work as before
  - Option names are checked against one list: names that differ only in upper and lower case are corrected (e.g. 'fontsize' -> 'FontSize'), 'Rotation' is accepted for `fontangle` (as documented; before, it was ignored), and unknown names give a warning once per session (identifier `mnt:netgraph:unknownOption`)
  - Colours: with matlab-utils from November 2026 on, the default colormap `rb_colors` goes from blue (negative values) to red (positive values)

## 1.3.0 (2026-10-31)

First version with a release tag (the previous version number was 1.2).

- `mnt_version` now gives three-part version numbers; `mnt_dependencies` checks required versions
- Bug fixes, SBML export, tests and README
- `mnt-extensions` renamed to `mnt-control-analysis`; `thermo_feasible_flux` added
- CPLEX replaced by `lp_solve` and `qp_solve` (optionally Gurobi)
- New default prior table `mnt_prior.tsv` (function `mnt_prior_file`); other packages can pass their own table to `parameter_balancing_prior`
- Moved here: `parameter_balancing_prior` (from Parameter balancing), `network_derivative` (from matlab-utils), `display_periodic_quantities_on_network` (from optimal-enzyme-rhythms)
- `state_data_to_data` moved to Model balancing
- `mnt_dependencies(1)` checks only the required packages
- Conversion between SBtab and network models, omics data, SBML and COBRA moved here from SBtab for Matlab (folder `mnt_sbtab_conversion`), including their tests
- `sbtab_load_network_model_and_data` uses `mnt_prior_file` as the default prior table (before: the prior table of Model balancing)
- `biotable_combine_replicates` moved to Replicate regression
- `es_default_options` renamed to `stm_default_options`; `es_check_feasibility`, `es_check_flux_constraints`, `es_check_velocities_ma_ratios`, `es_mu_simple_guess` and `es_make_fluxes_stationary` renamed to `mnt_...`; the structs `es_options` and `es_constraints` are now `stm_options` and `stm_constraints`
- Bug fixes:
  - `internal_elasticities` and `convenience2vector`: reactions in which a metabolite is both a reactant and a regulator, with equal values in `N` and `regulation_matrix`, were left out (because of a mistake in operator precedence), so `internal_elasticities` gave wrong results for convenience kinetics
  - `network_velocities` now respects the argument `indices` also for modular and mass-action kinetics; before, it returned all reactions, so `internal_elasticities` stopped with an error for these kinetics
  - `network_efmtool` no longer uses a fixed folder name and restores the current folder also after errors
  - `optimise_delta_mu_given_signs` and `optimise_mu_given_signs` work again (syntax error, missing helper functions)
  - `network_build_from_sum_formulae`: mistake in the check of the size of the regulation matrix
- `netgraph_move2` removed (it contained a syntax error and could not be called); `netgraph_move` remains
- The second copy of `EBA_orth` (in `eba_beard`) removed; the remaining copy is the one that was used before, with the usual Matlab path order
- `nanmean`, `nanstd`, `nanmax`, `nanmin`, `isstr` and `findstr` replaced by the recommended functions
- Functions that were not part of the toolbox moved here, so these features no longer depend on unpublished code: kinetics of type `standard` (`set_standard_kinetics` and related functions) and helper functions of `stochastic_simulations`, `display_MC_and_analytic`, `display_histograms` and `osc_resp_plot_timecourses`
- Kinetics of type `embedded_kinetic_models` give a clear error if `embedded_kinetic_velocities` is missing
- `ms_import_kinetic` deprecated; it gives a clear error if `ncsp_kinetic_data_load_table` is missing
- HTML documentation removed from the repository (it can be regenerated with m2html)
- Help texts added to functions that had none; `CITATION.cff` and policy for removing functions added

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
