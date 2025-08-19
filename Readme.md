
# habits1: C++ code for evolutionary simulation of habit forming and breaking


## Overview

This repository contains C++ code and example data.
The executable programs `EvoLearn`, built from this code, will run evolutionary simulations of a population of several groups of individuals that explore and form and break habits.
The program was used to produce results for the paper "Evolution of behavioural flexibility and the forming and breaking of habits" by Olof Leimar, Sasha R. X. Dall, Peter Hammerstein, Alasdair I. Houston, Bram Kuijper, and John M. McNamara.


## System requirements

The program has been compiled and run on a Linux server with Ubuntu 24.04 LTS.
The C++ compiler was g++ version 13.3.0, provided by Ubuntu, with compiler flags for c++20, and `cmake` (<https://cmake.org/>) was used to build the program.
It can be run multithreaded using OpenMP, which speeds up execution times.
Most likely the instructions below will work for many Linux distributions.
The program has also been compiled and run on macOS, using the Apple supplied Clang version of g++.

The programs read input parameters from TOML files (<https://github.com/toml-lang/toml>), using the open source `toml.hpp` header file (<https://github.com/ToruNiina/toml11/blob/main/single_include/toml.hpp>), which is included in this repository.

The programs also read compound stimulus data and gamete data from CSV files, using the open source `rapidcsv.h` header file (<https://github.com/d99kris/rapidcsv/blob/master/src/rapidcsv.h>), which is included in this repository.

The programs store evolving populations in HDF5 files (<https://www.hdfgroup.org/>), which is an open source binary file format.
The program uses the open source HighFive library (<https://github.com/BlueBrain/HighFive>) to read and write to such files.

These pieces of software need to be installed in order for `cmake` to successfully build the program.


## Installation guide

Install the repository from Github to a local computer.
There is a top directory `habits1`, containing this Readme.md file and a license file, and two subdirectories, OE and TS (Optimistic exploration and Thompson sampling), each with source code and a subdirectory `Data` where input data and data files containing simulated populations are kept (as well as subdirectories `build` used by `cmake` for files generated during building, including the executables `EvoLearn`.


## Building the programs

The CMake build system is used.
If it does not exist, create a build subdirectory in the project folder, OE or TS, (`mkdir build`) and make it the current directory (`cd build`).
If desired, for a build from scratch, delete any previous content (`rm -rf *`).
Run CMake from the build directory. For a release build:
```
cmake -D CMAKE_BUILD_TYPE=Release ../
```
and for a debug build replace Release with Debug.
If this succeeds, i.e. if the `CMakeLists.txt` file in the project folder is processed without problems, build the program:
```
cmake --build . --config Release
```
This should produce an executable in the `build` subdirectory.


## Running

Make the Data directory current (in either OE or TS).
Assuming that the executable is called `EvoLearn` and with an input file called `Run02h.toml`, producing detailed data for Figure 2b in the paper, run the program as
```
../build/EvoLearn Run02h_inp.toml
```
For evolutionary individual-based simulations, used for instance to determine evolutionary equilibria, the R script file `Run02_run.R`, can be used as follows
```
Rscript Run02_run.R
```
where `Rscript` is the app for running R scripts. You need to have `R` installed for this to work. Also, for this to work you need a starting population in a file `Run02_pop.h5`, which can be created by running
```
../build/EvoLearn Run02_inp.toml
```
with the input file `Run02_inp.toml` modified to contain `read_from_file = false`.


## Description of the evolutionary simulations

There is an input file, for instance `Run02_inp.toml`, for a case, which typically simulates 1000 years (generations), inputting the population from, e.g., the HDF5 file `Run02_pop.h5` and outputting to the same file.
Without an existing `Run02_pop.h5` data file, the program can start by constructing individuals with genotypes from the allelic values given in a gamete CSV file `Run02_gam.csv`in the input file.
To make this happen, use `read_from_file = false` in the input file.
Once you have a HDF5 file with a simulated population, there is an R script, e.g. `Run02_run.R`, which repeats runs a number of times, for instance 100, and for each run computes statistics on the evolving traits and adds a row to a CSV data file, e.g. `Run02_data.tsv`.

Each case is first run for many generations. When a seeming evolutionary equilibrium has been reached, 100 runs are kept in the summary data file, for instance `Run02_data.csv`. The equilibrium average allelic values are put into a CSV file, e.g. `Run02h_gam.h5`, which is then used for a simulation of detailed data and output to a file, e.g. `Run02_hist.h5` (such simulations are used for figures, e.g. for figure 2 in the paper).


## License

The `EvoLearn` program runs evolutionary simulations of habit forming and breaking.

Copyright (C) 2025  Olof Leimar

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
