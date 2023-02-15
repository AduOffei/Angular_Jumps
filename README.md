# Angular_Jumps
This repository contains the codes for characterizing the angular jumps in a molecular dynamics simulation of water


This file contains instructions for characterizing the angular jumps in a molecular dynamics simulation of water.

- The dataset file is called "prd_space.pdb"

- Software dependencies: matlab (tested in version 9.1), gfortran (tested in version f90) 

- There are no requiered non-standard hardware requirements

Instructions:

The procedure is organized in 4 steps, each one with a corresponding analysis tool.

Step 1: Extraction of the hh and dipole vectors 

code: dphh.f90 (Fortran90)

compilation: gfortran -o dphh.exe  dphh.f90

execution: "./dphh.exe"
Input: a trajectory file in .pdb format.

Output: the dipole and hhvectors of every molecule along the simulation time. 
The name of the outputs is DP_space.dat and HH_space.dat

Parameters (internal variables in dphh.f90)
 numframes=number of frames
 L = length of box
 kg = number of molecules

Expected runtime: ~2 minutes in a 2.3 GHz Intel Core i5.

Step 2: Angular jump detection protocol 

code: aj_detection.m (matlab)

input: - dp_space.dat or hh_space.dat
       - trajectory length
       - DP vector components =3x Number of water molecules 

output: - gd_angle_vs_time_dp_1.dat -> it has 5 columns: (i) the swing magnitude, (ii) swing duration, (iii) end time of swing, (iv) number of molecule

Parameters (internal variables in aj_detection.m)
t_step = inverse of your sampling frequency (time delay between snapshots in your trajectory file)
num_mol = number of molecules
optional
f_c = cutoff frequency (for low pass filter)
 

Step 3: Local topology extraction 

code: Defects.f90 (Fortran90)

input: the trajectory file in .pdb format

output: "defects_space.dat" contains the local topology of every molecule, i.e. the number of donors and acceptors.

parameters(internal variables to the code):
kg= number of molecules
numframes=number of frames
L=size of box

%Hydrogen bond definition
ang_ct= angular cutoff
dist_ct=Distance= distance cutoff


Expected runtime: ~7 minutes in a 2.3 GHz Intel Core i5.


Step 3: Nearest hydrogen bonded neighbors neighbors 

code: neighbors.f90 (Fortran90)

input: a pdb trajectory file 

outputs: "neigh_space.dat" the list of the neighbors for every water molecule (maximum donors =3, maximum acceptors = 3)

parameters(internal variables to the code. To be modified in the f90 file):
kg= number of molecules
numframes=number of frames
L=size of box

%Hydrogen bond definition
ang_ct= angular cutoff
dist_ct=Distance= distance cutoff

Expected runtime: ~1.5 minutes in a 2.3 GHz Intel Core i5.

Step 4: Nearest hydrogen bonded neighbors

code: nn.m (matlab)

Input: a trajectory in .pdb format, with frames 1ps apart from each other, and the output file of step 2 "gd_angle_vs_time_xx_1.dat"

output: the plots of the 1st 2nd 3rd 4th nearest neighbor distributions.


