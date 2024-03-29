# EPEC 

## Description

    Suite of programs to simulate multi-harmonic magnetic island dynamics in 
    presence of resonant magnetic perturbation  in time-varying toroidal  
    tokamak equilibrium.

## Contents

 ### Makefile
    GNU makefile for entire package
 
 ### /DOCUMENTATION
    - EPEC.tex: Latex description of compelete FLUX/NEOCLASSICAL/PHASE algorithm
 
 ### /SCRIPTS
    Useful shell scripts

 ### /FLUX 
    Program to read gFile(s) and write equilibrium data needed by NEOCLASSICAL 
    and PHASE to fFile

 ### /NEOCLASSICAL 
    Program to read FLUX data from fFile(s), profile data from pFile(s) and cFiles(s), 
	and write neoclassical data needed by PHASE to nFile

 ### /PHASE 
    Program to read FLUX data from fFile(s), NEOCLASSICAL data from nFile(s), 
	GPEC data from uFile(s), mFiles(s), and lFile(s), and perform island dynamics 
	simulation in fixed equilibrium. Final state of plasma saved in sFile.

 ### /RESCALE 
    Program to rescale equilibrium gFile to modify q_95 by modifying toroidal 
    plasma current while keeping B_toroidal the same

