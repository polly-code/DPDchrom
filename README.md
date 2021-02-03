# DPDchrom - dissipative particle dynamics for chromosome simulations

## General description

DPDchrom is developed to reconstruct 3D chromatin conformation using single cell Hi-C contact map. This method is based on DPD [1]. Previous version with manual preparation published [2] and available [here](https://github.com/polly-code/DPD_withRemovingBonds). It means that thermostat is embedded into a motion calculation. Solvent is taken into account explicitly. Due to the soft repultion, one can use large integration time step, dt=0.04. This software takes contact matrix as an input and produces chromatin conformation as an output.

## Installation

To install DPDchrom, put in the same folder makefile and DPDchrom.f90, then compile the fortran code with any fortran compiler, simply type 
> make

DPDchrom is ready to run.

## How to use

In the test folder you can find example of the input **37420_*.csv** [3]. The format of the file is: 

> chr1,chr2,start1,end1,start2,end2,count

> 1,1,3200000,3400000,3200000,3400000,1

To perform calculations execute DPDchrom in the command line passing two arguments: "/path/to/file" and resolution in bp "100000". Example of the command to reconstruct conformation using file **myfile** at resolution **100 Kb**:

> ./DPDchrom /home/work/myfile 100000

## Postprocessing

After simulation you will have your structure in unknown format. You can use script **rst2mol2.py** from py_script folder to convert it into **mol2** format. Moreover you can also remove periodic boundary conditions properly.

To compare your structures with each other, please, use script **imj_acc.py** as template. It allows you to calculate accuracy (Modified Jaccard Index).

## Spin-off

In the folder cmd you can find input for CMD calculations for the LAMMPS.

## References

1. Groot, Robert D., and Patrick B. Warren. "Dissipative particle dynamics: Bridging the gap between atomistic and mesoscopic simulation." The Journal of chemical physics 107.11 (1997): 4423-4435.
2. Ulianov, Sergey V., et al. "Order and stochasticity in the folding of individual Drosophila genomes." Nature Communications 12.1 (2021): 1-17.
3. Flyamer, Ilya M., et al. "Single-nucleus Hi-C reveals unique chromatin reorganization at oocyte-to-zygote transition." Nature 544.7648 (2017): 110-114.
