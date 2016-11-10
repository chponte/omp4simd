### Evaluation of OpenMP 4 SIMD Directives on Xeon Phi Coprocessors

#### Introduction

This repository contains the source code for all the benchmarks used in the study as well as the scripts utilized for the results gathering.

#### Compilation

After downloading the appropiate version of the codes, a configuration needs to be used. Replace the `config.make` configuration with one of the configuration templates from the configs folder, and adapt it to your needs. Lastly, run make to generate the binary files.

```bash
$ cd omp4simd
$ cp configs/config.autovec config.make
$ vi config.make
$ make
```

#### Usage

You can either copy the generated  binaries onto the Intel Xeon Phi and manually execute them or use the scripts included in the src/scripts folder to automatice the result gathering.

In order to use the scripts, `screen` needs to be installed on your system. Before running any of them, it is necessary to adapt the definitions included in each file to your system. After all the changes are done, simply run `src/scripts/runmic.sh`. The usage is as follows: `./runmic.sh <address> <path to program folder> <path to results folder> <repetitions per program> [label]`. You can check wether the execution is still in progress by running `screen -ls`. Results are obtained in `csv` format, ready to be imported into a spreadsheet.

```bash
$ cd omp4simd
$ src/scripts/runmic.sh 0.0.0.0 bin/ results/ 10 "default execution"
$ screen -ls
```

##### Definitions

runmic.sh
```bash
# MIC libraries path
LIB_PATH=path_to_the_intel_mic_libraries
# Task script
TASK=path_to_task.sh_script
```

task.sh
```bash
# Scratch directory
SCRATCH_DIR=path_to_scratch_directory
# Thread number
THREAD_NUMBER=number_of_threads_during_execution
```

#### Help

For questions, suggestions or problems please use the Issues tab here.

#### Licensing

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.