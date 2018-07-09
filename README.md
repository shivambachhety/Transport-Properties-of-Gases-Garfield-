# Transport-Properties-of-Gases-using Garfield++
To simulate the transport properties of electrons in gas mixtures and visualize Paschen curves
## Introduction
I have used a gas mixture of 90% argon and 10% iso-butane at a pressure of 3 atm and room temperature to study transport parameters (drift velocity, diffusion coefficients, Townsend coefficient, and attachment coefficient) as a function of the electric field E (and, in general, also the magnetic field B as well as the angle between E and B) and using SRIM (Stopping and Range of Ions in Matter) for simulating the energy loss of ions in matter and obtaining plots for stopping powers, range and straggling parameters.

## Prerequisites
* Garfield++: An object-oriented toolkit for the detailed simulation of particle detectors that use a gas mixture or a semiconductor material as sensitive medium.
* ROOT (v 5.34): A modular scientific software framework. It provides all the functionalities needed to deal with big data processing, statistical analysis, visualisation and storage
* C++ compiler that supports C++11 and a Fortran compiler

## Simulation and Plots
I have carried out simulation of Argon- isobutane (90%-10%) gas mixture using 20 electric field points between 100 V/cm and 100 kV/cm with logarithmic spacing and using Magboltz to specify the number of collisions (in multiples of 107). We have also obtained values for different parameters using monte-carlo simulation, Paschen curve and the energy plots using SRIM.
 
## Conclusion
Therefore, I carried out the simulations for a different mixture of gases under the electric and magnetic field and the corresponding gas table is saved in ar_90_ic4h10_10.gas. The gasfile.C can be edited for different values of electric fields

## Results
The results of gas tables (transport parameters) as shown in Fig. below and calculating electric fields can be used in modelling a drift tube.

 

## References
1.	Garfield++ – simulation of tracking detectors. https://garfieldpp.web.cern.ch/garfieldpp/
2.	ROOT- https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html
3.	Mobility and diffusion of ions in gases- McDaniel, E.W.; Mason, E.A.
4.	Measurements and simulations of drift gas properties by Lukas Koch, http://www.institut3b.physik.rwthachen.de/global/show_document.asp?id=aaaaaaaaaakafeg
5.	https://garfieldpp.web.cern.ch/garfieldpp/examples/gem/GemTransport.pdf
