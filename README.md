# Chemical Equilibrium of Glycerol Conformers

This repository contains Fortran90 codes to calculate the partition function and equilibrium constant between two conformers of glycerol (A and B as shown in image "immagine_conformeri.JPG")), using "ab initio" data  including vibrational frequencies and zero-point energy in the gas phase.
The data are extracted from the paper:

------ J. Phys. Chem. A 2000, 104, 47, 11220–11222, Riccardo Chelli et al., [DOI: 10.1021/jp002677e](https://doi.org/10.1021/jp002677e) ------

The calculation is also focused on the temperature dependence of the equilibrium constant.

===== CODE STRUCTURE =====

       --M-A-I-N--        
      
  - `funz_partizione.f90`: Main program for calculating the partition function at T = 298.15 K.
  - `funz_partizione_temperature.f90`: Main program for partition function over a range of temperatures (40 K to 420 K).
 
    --M-O-D-U-L-E-S--
                                        
  - `Vibrazionale.f90` and `mom_inerzia.f90`: Modules for calculations at T = 298.15 K.
  - `Vibrazionale_temperature.f90` and `mom_inerzia_temperature.f90`: Modules for range of temperatures.
    
====== Additional Scripts =======

  - A basic MATLAB script is provided to fit the equilibrium constant vs. temperature data and visualize the results(As shown in Image kp_T)

 ==== Results ====            
 
The Code provide that conformer B is more stable than A, at 298.15 K and generally all over the range.
The analysis reveals an exponential decrease in the equilibrium constant with increasing temperature. This behavior is visualized using the provided MATLAB script.

=== How to Compile ====
1. --- Compile Modules ---
   gfortran -c Vibrazionale.f90 mom_inerzia.f90
2. --- Compile Main ---
   gfortran -o main funz_partizione.f90 Vibrazionale.o mom_inerzia.o
   /main.exe
The same for the code over the temperature
