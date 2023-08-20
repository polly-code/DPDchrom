# DPDchrom - dissipative particle dynamics for chromosome simulations

## General description

DPDchrom is developed to reconstruct 3D chromatin conformation using a single-cell Hi-C contact map. This method is based on DPD [[1](#references)]. Previous version with manual preparation published [[2](#references)] and available [here](https://github.com/polly-code/DPD_withRemovingBonds). It means that the thermostat is embedded into a motion calculation. A solvent is taken into account explicitly. Due to the soft repulsion, one can use a large integration time step, dt=0.04. This software takes the contact matrix as an input and produces chromatin conformation as an output.

## Installation

To install DPDchrom, put in the same folder `makefile` and `DPDchrom.f90`, then compile the Fortran code with any Fortran compiler. 
> make

DPDchrom is ready to run.

## How to use

In the test folder there is an example of input `37420_*.csv` [[3](#references)]. The format of the file is: 

`chr1,chr2,start1,end1,start2,end2,count`

`1,1,3200000,3400000,3200000,3400000,1`

To perform calculations execute DPDchrom in the command line passing **3** arguments, _/path/to/file_, resolution in bp, and chain length. Path to file and resolution are self-descriptive. If you want to reconstruct a single chain with a fixed length, then pass the number of beads. When you don't want to define chain length or have several chains, please, pass **0**. Example of the command to reconstruct conformation using file `myfile` at resolution _100 Kb_:

> ./DPDchrom /home/work/myfile 100000 0

## Postprocessing

The format of the output file with conformation is custom. In order to convert it to _mol2_ format, there is a script `rst2mol2.py` in _py_script_ folder. Moreover, there is an opportunity to remove periodic boundary conditions properly.

To compare structures with each other, there is a script `imj_acc.py` provided as a template to calculate similarity (reconstruction accuracy) as the Modified Jaccard Index.

## Spin-off

In the folder _cmd_ you can find input files for CMD calculations for [LAMMPS](https://github.com/lammps/lammps).

## References

1. Groot, Robert D., and Patrick B. Warren. "Dissipative particle dynamics: Bridging the gap between atomistic and mesoscopic simulation." The Journal of chemical physics 107.11 (1997): 4423-4435.
2. Ulianov, Sergey V., et al. "Order and stochasticity in the folding of individual Drosophila genomes." Nature Communications 12.1 (2021): 1-17.
3. Flyamer, Ilya M., et al. "Single-nucleus Hi-C reveals unique chromatin reorganization at oocyte-to-zygote transition." Nature 544.7648 (2017): 110-114.
