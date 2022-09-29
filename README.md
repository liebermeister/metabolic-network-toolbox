Metabolic Networks Toolbox for Matlab
=====================================

The Metabolic Networks Toolbox contains functions for constructing, simulating, and analysing computational models of metabolism.

## Dependencies:

  o Matlab utility functions (https://github.com/liebermeister/matlab-utils)

  o SBML toolbox    (http://sbml.org/Software/SBMLToolbox)

  o SBtab toolbox  (https://github.com/liebermeister/sbtab-matlab)

  o Tensor toolbox (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)

  o efmtool        (http://www.csb.ethz.ch/tools/efmtool)

For reading and writing cobra files, the matlab cobra toolbox is also required. Please make sure that these matlab packages are installed in your system and that all these directories and subdirectories are included in your matlab path.

Some functions employ linear or quadratic optimisation. If CPLEX is installed in your system, version 12.10 is required. If CPLEX is not installed, native matlab functions for optimisation are used instead.  

## Installation
Please see the [installation instructions](INSTALLATION).

HTML documentation is available in the subdirectory doc

## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.
