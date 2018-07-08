# tish
![version][version-image]  
Waveform computation by *Direct solution method* (DSM)

This software is based on the software which can be downloaded from [University of Tokyo][dsm]

**Welcome to the DSM world.**

This is a package for SH computation 
for transversely isotropic spherically symmetric media 
in a frequency domain.
You need FORTRAN90 to run this program.


CONFIGURATION:  
```	bash
% make
```
USAGE:  
``` bash
% mpi-tish < (parameter filename)
```	

This *mpi-tish* can also run by mpi.

[dsm]: http://www-solid.eps.s.u-tokyo.ac.jp/~dsm/software/software.htm
[version-image]:https://img.shields.io/badge/version-0.90.0-yellow.svg
