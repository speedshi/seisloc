# seisloc

Seismic location using waveform migration based methods,
including the Multichannel Coherency Migration method of Shi et al. (2018a, 2018b).

The codes can be freely used for NonCommercial activities, such as research and teaching.

## Installation

To install the `migrationloc` programs, you need:

- A Fortran compiler
- An MPI installation for the MPI version of the code

Build the software by changing to the root directory of the repo, and typing:

    make -C src

By default this will build the serial, OpenMP and MPI versions.  Other `make` targets:

- `serial` : Only build the serial version
- `mpi`: Only build the MPI version
- `openmp`: Only build the OpenMP version

### Compiler and compilation flags

By default, we use `gfortran` to compile the programs and `mpif90` to link the MPI libraries.

To use the Intel compiler, set `FC=ifort` and `MOD=module`.


## References

More details about the theory and applications of MCM can be found in:

1. Shi, P., Angus, D., Rost, S., Nowacki, A. and Yuan, S., 2018a. Automated seismic waveform location using Multichannel Coherency Migration (MCM)–I. Theory. _Geophysical Journal International_, 216, 1842&ndash;1866.
doi:[10.1093/gji/ggy132](https://doi.org/10.1093/gji/ggy132)

2. Shi, P., Nowacki, A., Rost, S., Angus, D., 2018b. Automated seismic waveform location using Multichannel Coherency Migration (MCM)–II. Application to induced and volcano-tectonic seismicity. _Geophysical Journal International_, 216, 1608&ndash;1632.
doi:[10.1093/gji/ggy507](https://doi.org/10.1093/gji/ggy507)

Please cite these publications if you use the MCM in your work.


## `src/envhlbfft.F90`

The FFT implementation used in `migrationloc` is distributed under the terms of the [GNU General
Public License](https://www.gnu.org/licenses/old-licenses/gpl-1.0.html).

```
License:

  Licensed under the GNU General Public License (GPL).
  Copyright (C) 1995-2004, Scientific Computing Division,
  University Corporation for Atmospheric Research

Author:

  Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
  FORTRAN90 version by John Burkardt.

Reference:

  Paul Swarztrauber,
  Vectorizing the Fast Fourier Transforms,
  in Parallel Computations,
  edited by G. Rodrigue,
  Academic Press, 1982.

  Paul Swarztrauber,
  Fast Fourier Transform Algorithms for Vector Computers,
  Parallel Computing, pages 45-63, 1984.
```