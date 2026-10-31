# Conventions for this toolbox

These conventions hold for new code and for files that are changed anyway. Existing files
do not need to be rewritten just to follow them.

## Compatibility

Other packages and scripts depend on this toolbox, so function names and fields of data
structures are public. Folders are not: code is found through the Matlab path, so files can be moved.

- Before renaming or removing a function, search for its callers in other packages.
- Follow the rule in `README.md`: an old name remains as a wrapper with a warning for at least
  one MINOR version, and is removed only in the next MAJOR version.
- List each change in `CHANGELOG.md`.
- A new version number must be set in `mnt_version`, `mnt/Contents.m` and `CITATION.cff`
  (`tests/test_version.m` compares the first two).

## Function names

| Topic | Convention | Example |
|---|---|---|
| Word order | object first, then verb | `state_data_save`, not `save_state_data` |
| Conversions | `x_to_y`, not `x2y` | `network_to_sbtab` |
| Matlab data structures in files | `load` and `save` | `kinetic_data_load` |
| External file formats | `import` and `export` | `network_sbml_export` |
| Language | British spelling | `optimise`, `analyse` |

## Function arguments

- A function has few positional arguments. Further settings go into a struct called `options`.
- Default values of options are defined in the function, or in `*_default_options` if users need to
  see or edit them, and are combined with `join_struct`.
- Default values of arguments use `eval(default('name','value'))`, like in most functions of the toolbox.
- Metabolite vectors have the order of `network.metabolites`, and reaction vectors the order of
  `network.actions`. Vectors are column vectors.

## Help texts

Each function starts with a help text in this form (the first line is used by `lookfor` and `help`):

```matlab
function [v, value] = my_function(network, s, options)

% MY_FUNCTION One-line description, without a full stop at the end
%
% [v, value] = my_function(network, s, options)
%
% Longer description, if needed.
%
% Inputs:
%   network   network structure (see 'help metabolic_network')
%   s         metabolite concentrations (column vector)
%   options   (optional) struct with fields:
%     verbose   (default 0) flag: display progress
%
% Outputs:
%   v         reaction rates (column vector)
%   value     ...
%
% Example:
%   ...
%
% See also: other_function, another_function
```

Examples must work when copied to the Matlab command window.

## Messages

- Functions print nothing unless `options.verbose` (or a `verbose` argument) is set.
- Error messages say what went wrong and what to do.
- Errors and warnings have an identifier: `error('mnt:missingField', 'The network has no field ...')`.

## Figures

Graphics functions draw into the current axes and do not choose figure numbers themselves.

## Code layout

- Indentation by 2 spaces, no tabs, no spaces at the end of lines
- `&&` and `||` in conditions; never `+` or `*` to combine conditions, because of operator precedence
- Bug fixes come with a test in the folder `tests`
