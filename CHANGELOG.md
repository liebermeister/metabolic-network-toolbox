# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.3.0 (not released yet)

First version with a release tag (the previous version number was 1.2).

- `mnt_version` now gives three-part version numbers; `mnt_dependencies` checks required versions
- Bug fixes, SBML export, tests and README
- `mnt-extensions` renamed to `mnt-control-analysis`; `thermo_feasible_flux` added
- CPLEX replaced by `lp_solve` and `qp_solve` (optionally Gurobi)
- New default prior table `mnt_prior.tsv` (function `mnt_prior_file`); other packages can pass their own table to `parameter_balancing_prior`
- Moved here: `parameter_balancing_prior` (from Parameter balancing), `network_derivative` (from matlab-utils), `display_periodic_quantities_on_network` (from optimal-enzyme-rhythms)
- `state_data_to_data` moved to Model balancing

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
