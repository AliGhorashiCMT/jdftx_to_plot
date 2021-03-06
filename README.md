[![Build status][ci-status-img]][ci-status-url]


# jdftx_to_plot
Plot functionalities and dielectric properties from JDFTX. 

## Installation

The package is written in Julia and may be added via 
```julia
pkg> add https://github.com/AliGhorashiCMT/jdftx_to_plot
```
All functionality may then be accessed by
```julia
julia> using jdftx_to_plot
```
Current capabilities are: Density of states calculations, 2d and 3d polarization/dielectric function calculations, phonon and electronic band structure calculations, electron-phonon matrix element calculations, momentum matrix element calculations, JDFTX input file creation, surface plasmon loss calculations up to 2nd order in phonon-assisted processes. In addition, example files are included in the package.

Throughout the package, the following conventions are used: Distance is in angstrom, wavevectors are in inverse angstrom, energies are in eV, time is in seconds


[ci-status-url]: https://travis-ci.com/github/AliGhorashiCMT/jdftx_to_plot/builds/222201217
[ci-status-img]: https://travis-ci.com/AliGhorashiCMT/jdftx_to_plot.svg?branch=main
