# DPDchrom

## Installation

Compile the fortran code with any fortran compiler. In the repository, you can find the Makefile to compile the program. Then DPDchrom is ready to run.

## How to use

After compiling DPDchrom, pass file path and desired resolution in bp. The format of the file is: 

> chr1,chr2,start1,end1,start2,end2,count
> 1,1,3200000,3400000,3200000,3400000,1

Example of the command to reconstruct conformation using file **myfile** at resolution **100 Kb**:

> ./DPDchrom /home/work/myfile 100000
