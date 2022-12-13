# Diagrammatic Monte Carlo Approach to Angular Momentum in Quantum Many-Particle Systems

This repository contains the code developed for the paper G. Bighin, T. V. Tscherbul, and M. Lemeshko, Phys. Rev. Lett. **121**, 165301 (2018). We introduce a novel approach to describe angular momentum in a quantum many-body system via Diagrammatic Monte Carlo. Please refer to [the paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.165301), also available on the [arXiv](https://arxiv.org/abs/1803.07990), for more details about the physics of the system.

# Requirements

A modern C compiler (GCC or Clang) and the following libraries:

- [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)
- NCurses

The following libraries/codes are included in the repository:

- [Cubature](http://ab-initio.mit.edu/wiki/index.php/Cubature_(Multi-dimensional_integration)), for adaptive multidimensional integration.
- [inih](https://github.com/benhoyt/inih) for parsing .ini files.
- [libprogressbar](https://github.com/doches/progressbar) to display a nice progress bar.
- MurmurHash3, a non-cryptographic hash
- [gnuplot_i](https://github.com/mithodin/gnuplot_i), a gnuplot interface in ANSI C
- spline.c for spline interpolation

# Usage

Compile the code using the provided Makefile (`make`), eventually adapting it to point to the correct location of your C compiler and of your libraries. Prepare a .ini file with the physical details of the system you want to simulate, the `diagmc.ini` contains all the possible options, along with comments. Finally run the code (`./diagmc diagmc.ini`). Depending on the options specified, a progress bar might or might not be shown while the program is running.

# Brief description of the repo

The entry point is in `main.c`, whereas the main Monte Carlo loop is in `mc.c`, function `do_diagmc()`; it also contains the basic high-level update logic, see the functions `update_length()`, `update_add_phonon_line()`, `update_remove_phonon_line()`, `update_shuffle()`.

The lower level update logic is contained in `updates.c`, working directly on the `struct diagram_t` as defined in `diagrams.h`, also note the very important function `diagram_weight()` in `diagrams.c`. On the other hand `graphs.c` contains routines used for simplifying diagram summations using ideas from the diagrammatic theory of angular momentum, see for instance the last chapters in D.A. Varshalovich, A. N. Moskalev, V.K. Khersonskii, *"Quantum Theory of Angular Momentum"* (World Scientific).


The files `histograms.c` and `stat.c` contain simple routines used for statistics and histogram samplings.

The files `phonon.c`, `physics.c` and `selfenergies.c` implement the physical routines used in calculating the diagram weights.

The file `aux.c` contains various auxiliary routines, while `config.c` reads configuration files.

Finally, `debug.c` and `tests.c` implement some debugging and testing infrastructure.

The folder `slurm` contains scripts to run the code on a SLURM cluster. The folder `xtra` contains experimental, untested code.
