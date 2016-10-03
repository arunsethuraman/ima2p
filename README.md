# IMa2p
IMa2p is a parallel implementation of IMa2, using OpenMPI-C++ - a Bayesian MCMC based method for inferring population demography under the IM (Isolation with Migration) model.
Please refer to Sethuraman and Hey (2015) for details of implementation.

Citation:
Sethuraman, A, J Hey. 2015. IMa2p - Parallel MCMC and inference of ancient demography under the Isolation with Migration (IM) model. Molecular Ecology Resources
DOI: 10.1111/1755-0998.12437

Important - Bug fix - 8/4/2016

Fixed a bug which was causing the program to hang in the event that the user doesn't specify the heating parameters despite trying to run it in parallel.
Thanks to Jared Knoblauch for catching this!

Important - Bug fix - 10/19/2015

Fixed a bug in reading MCF files back into memory for re-starting M mode. Previously, the code was wrongly reading the MCF file prefixed as:
"newoutputfilename.mcf.processorid", whereas it should have been "oldoutputfilename.mcf.processorid". I have now fixed this. You should be able to specify the file name in command line as:

mpirun -np ... -r3 -f oldoutputfilename.out

Thanks to Louis Plough for bringing this to my attention!


Important - Bug fixes - 8/3/2015

Users that downloaded IMa2p before 7/20/2015, please note that there have been a couple of minor bug fixes since.

1) TI files - all genealogies sampled by IMa2p will be printed to a file with a .ti extension. The buggy version was printing files suffixed with .ti.0 instead, which in turn was leading to errors while calling the L mode.
This is a simple fix - you will have to rename the file to have a .ti extension instead, and then rerun IMa2p L mode.
To do this in Unix, type:

mv yourfile.ti.0 yourfile.ti

Thereon, you should be able to call IMa2p's L mode as before. Note that this has now been fixed in the new version, and the program will directly print .ti files instead.

2) I have now added the functionality to print swap rates between adjacent temperatures - note that in parallelized MC3, while chains remain stationary across processors, their temperatures move around throughout the run.
Hence swapping rates between adjancent temperatures of chains should give a better idea of mixing than swap rates between adjacent chains. 

Please report any bugs/crashes to asethuraman@csusm.edu

Compiling and Running IMa2p

By Arun Sethuraman and Jody Hey

Center for Computational Genetics and Genomics, Temple University

1.	Introduction
This document details how to compile and execute IMa2p, a parallel version of IMa2, implemented under an OpenMPI-C++ framework. 
IMa2p uses the same command set as IMa2, and users should refer to the documentation on using IMa2 and IM/IMa for further details.

2.	Compiling using autoconf
IMa2p has been written using the OpenMPI-C++ framework, and can be compiled using standard MPI flavors of the GNU compiler (including mpicc and mpicxx).  OpenMPI must be installed.  For details on installing OpenMPI, please see this page. Also, once installed, please make sure that the “bin” directory of OpenMPI (which contains all the compilers, and executables) is added to your machine’s path.
If you have already installed OpenMPI, to compile using autoconf (assuming you have some stable version of autoconf on your Unix machine), download the package, unzip/untar it, then cd into the folder. Once inside make sure that you have permissions to install using:

chmod +x *

Then, type:

./configure --with-mpi=yes --prefix=/path/to/install

This should create the necessary ‘make’ files to install the package. 
If you aren’t sure, the configuration script will determine if your machine contains OpenMPI definitions or not, and create ‘make’ files accordingly. To do this:

./configure --with-mpi=auto --prefix=/path/to/install

Alternately, IMa2p can also be compiled using a generic GNU C++ compiler for serial use (similar to the original IMa2 package). To do this, type:

./configure --with-mpi=no --prefix=/path/to/install

Note however, that if your machine does not have OpenMPI installed, or if the configuration script is unable to find it, it will compile a serial version of the program (similar to IMa2) with a regular GNU C++ compiler.
Once the ‘make’ files have been generated successfully, you should be able to install IMa2p by typing:

make

This should create an executable called “IMa2p” inside the src folder. Alternately, you can create a separate executable folder by typing:

make install

This should create a bin folder inside your installation path, which will contain the executable (“IMa2p”).

3.	Compilation from source
If you wish to compile from source, then you need to either delete or comment out line #24 from src/imamp.hpp (which says #include <config.h>).
Thereon, to compile from source, cd into the src folder, and type:

mpicxx -DMPI_ENABLED –o IMa2p *.cpp

This will create an executable called IMa2p in the src folder. To run the program in serial mode (single processor), run the executable directly with an appropriate command line (for details see the manual for IMa2).
To run it in parallel use mpirun or mpiexec, for example:

mpirun –np  <number of processors to use> IMa2p … <IMa2 command line options>

An IMa2p run must invoke at least 1 MCMC chain per run.  In IMa2 the –hn flag is used to indicated the total number of chains.  In IMa2p the –hn flag is used to indicate the number of chains per processor.   For example, to invoke Metropolis coupling using a total of 10 chains, distributed among 5 processors:

mpirun –np 5 IMa2p  -hn 2 <other IMa2 command line options> 

For instructions on compilation on a multi-core Windows PC using a Linux Emulator (eg. Cygwin), please refer to the user manual.

4.	Information in the Output files: Differences  from IMa2:
a)	Main output file:   In addition to information on swapping rates between chains, the output file summarizes the swapping rate between processors. 
b)	*.ti file – if you have set genealogies to be saved into separate *.ti files, these will be created for each processor, but the sampled genealogies will be saved only on the head node. This file will be named *.ti.
c)	*.mcf.<processor number> files – if you have set the MCMC states to be saved on each chain, this information is saved in files with the *.mcf extension. Each chain on each processor saves a different state file. All these files are important, if and when you wish to restart your ‘M’ or ‘L’ mode run in parallel.
d)	*.burntrend.out.<processor number> file – burn-in trend files are created at the end of burn-in runs. These contain update rates for genealogies, which are unique to each processor (since there are coupled chains running on each processor). However, the burn-trend is collated onto the head node, and can be seen in the *.burntrend.out.0 file.
e) If you wish to save lots of genealogies from multiple separate runs (M mode), and then run L mode computations on all these genealogies, you can do so by using the -r0 option, and specifying suffixes for all the saved ".ti" files using the -v option. Note that a separate L mode run will run in serial.

5.	Please report any issues with compilation/running IMa2p to Arun Sethuraman (asethuraman@csusm.edu).


