# godunov2d_distributed_SOU_ArtViscosity
Solve 2d Euler equations using Godunov FVM. Utilized second-order upwind scheme with minmod flux limiter (ala MUSCL approach). 
Flux splliting is based on AUSMup or ROE methods. Time marching is possible by explicit first-order Euler or TVD RK3 schemes.

Added: 
artificial viscosity method for curing the shock instabilities (Rodionov, C&F2021; Rodionov);
Gradient, divergence and Laplace opertaors; 
van Leer, van Albada and modified minmod limiters. 

The test case is A Mach 3 Wind Tunnel With a Step.
For more than fifty years supersonic Mach 3 flow over forward step is a canonical benchmark to test different computation schemes. 
It was introduced for the first time by Emery, then was represented by Van Leer, and was popularized most widely by Woodward and Colella (1984). 
Special interest is general capability of numerical methods to reproduce complex aerodynamic phenomena, 
such as transient interactions of shock and rarefaction waves as well as Mach disks arising in the process of irregular interactions
of waves with each other and with the channel walls. It is well known that the corner of the step is a singularity point, 
which can produce significant artificial numerical errors depend on applied boundary conditions and numerical schemes used.

Figure below demonstrates solution at time, t= 4.0/sqrt(gamma),  obtained on the baseline 80x20 trinagular grid(s) with the smoothed corner:
contours of density (top) and artifical viscosity (bottom),  respectevily. 

![alt text](https://github.com/dimaZloy/godunov2d_distributed_SOU_ArtViscosity/blob/main/results2.png)
