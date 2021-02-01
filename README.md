# DPDchrom - dissipative particle dynamics for chromosome simulations

## General description

DPDchrom is developed to reconstruct 3D chromatin conformation using single cell Hi-C contact map. This method is based on DPD [1]. That means that thermostat is embedded into the calculation of forces. Solvent is used explicitly. Due to the soft repultion, one can use large integration time step. This software takes contact matrix as input and creates chromatin conformations as an output.

## Installation

To install DPDchrom, put in the same folder makefile and DPDchrom.f90, then compile the fortran code with any fortran compiler, simply type 
> make
DPDchrom is ready to run.

## How to use

After compiling DPDchrom, pass file path and desired resolution in bp. The format of the file is: 

> chr1,chr2,start1,end1,start2,end2,count
> 1,1,3200000,3400000,3200000,3400000,1

Example of the command to reconstruct conformation using file **myfile** at resolution **100 Kb**:

> ./DPDchrom /home/work/myfile 100000

## Postprocessing

After simulation you will have your structure in unknown format. You can use script **rst2mol2.py** from py_script folder to convert it into **mol2** format. Moreover you can also remove periodic boundary conditions properly.

To compare your structures with each other, please, use script **imj_acc.py** as template. It allows you to calculate accuracy (Modified Jaccard Index).

## Spin-off

In the folder cmd you can find input for CMD calculations for lammps engine.
