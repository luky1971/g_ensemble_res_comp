### g_ensemble_res_comp

g_ensemble_res_comp evaluates the difference between the residues of two conformational
ensembles, R and R'. Quanitification is in terms of a true metric that
satisfies the conditions set forth by the zeroth law of thermodynamics. The
quantification metric eta=1-|Overlap|=|R'|-|Overlap|=DeltaR is normalized,
that is, 0<=eta<1, and takes up a value closer to unity as the difference
between the ensembles increases.

The two ensembles are provided as two trajectory files specified
by the `-f1` and `-f2` options (only pdb files are supported).
We recommend that frame numbers in the trajectory files are in the range
2500-5000. While the speed of the algorithm decreases with increase in
ensemble size, the numerical accuracy of the calculation reduces with
decrease in ensemble size, and a small number of frames may not provide a good
representation of the ensemble.

By default, differences (eta) are estimated for all residues.
Overlaps are estimated by training a support vector
machine in a pre-defined Hilbert space specified by the width of the RDF
Kernel (gamma=0.4) and the maximum value that can be taken up by the
Lagrange multiplier (C=100.0).

By default, g_ensemble_res_comp is parallelized with OpenMP (see installation instructions below). The default behavior is to use the maximum number of cores available.

### INSTALLATION

The following instructions are for unix-based operating systems such as OSX and Linux. Windows is not currently supported.

1. Install Gromacs version 4.5.x or later from http://www.gromacs.org.

2. `git clone` or otherwise obtain and `cd` to the 'g_ensemble_res_comp' repository.

3. Run `sudo make install` with the necessary arguments for your environment (see below).

If you do not have Gromacs version 5.x installed, you will need to set the makefile's `VGRO` variable to the root Gromacs version number. If Gromacs is installed in a non-default directory (ie not in /usr/local/gromacs) then you will have to set the `GROMACS` variable to the Gromacs installation directory that contains the 'include' and 'lib' folders. For example, if you are running Gromacs 4.5.3 installed in /home/user/tools/gromacs-4.5.3, then you would run the following command:

``` bash
$ sudo make install VGRO=4 GROMACS=/home/user/tools/gromacs-4.5.3
```

If you must run `make install` without sudo privileges, you will need to set the `INSTALL` variable to a path that you can write to.
The default install path is /usr/local/bin. Depending on your system and chosen installation directory, you may have to add g_ensemble_res_comp to your PATH.

You will also need to make sure that `CC` and `CXX` are set to the C compiler and C++ compiler commands, respectively, that were used to compile your installation of Gromacs. For example, if your environment's default compiler is clang, but you used gcc to compile Gromacs, you would run the following:

``` bash
$ sudo make install CC=gcc CXX=g++
```

If you want to build without OpenMP, set `PARALLEL=0`. You can also add compilation flags by setting `CFLAGS`, and linker flags/libraries by setting `LIBS`. For example, if you set `LIBS=-static` to statically link g_ensemble_res_comp's dependencies, you can then run the same binary in a different environment without the same C runtime or Gromacs library present.

### USAGE

``` bash
$ g_ensemble_res_comp -f1 first_file.pdb -f2 second_file.pdb
```

if you run just g_ensemble_res_comp the default PDZ2_frag_apo.pdb and PDZ2_frag_bound.pdb will run.
