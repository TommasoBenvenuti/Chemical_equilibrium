# Chemical_equilibrium
In this Fortran90 code I calculated the partition function for two conformers of glycerol, starting from "ab initio" data, such as vibrational frequencies and zero-point-energy.
This data are extracted from : J. Phys. Chem. A 2000, 104, 47, 11220–11222, Riccardo Chelli et al. , https://doi.org/10.1021/jp002677e
I also tried to visualize the effect of temperature on the equilibrium. An exponential decrease of the equilibrium costant was expected, but in this way I managed to demonstrate this. 

I wrote a basic matlab script to fit this data, using fitting libraries. 

"Vibrazionale.f90", "mom_inerzia.f90" are modules to be compiled before the main "funz_partizione.f90". They are referred to the calculation at temperature T= 298.15 K.

"Vibrazionale_temperature.f90", "mom_inerzia_temperature.f90" and "funz_partizione_temperature.f90" are the codes referred to the calculation over the temperatures. 
