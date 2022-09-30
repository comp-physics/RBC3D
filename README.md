# Spectral Boundary Integral Solver for Electrohydrodynamic Flows in Viscous Drops
## Authors: M. Firouznia, S. H. Bryngelson, D. Saintillan

EHD_Drop_3D is an electrohydrodynamic solver for viscous drops. It solves the boundary integral forms of the Laplace's equations and Stokes equations via an algorithm tailored for drop-scale simulations:

* Spectrally-accurate spherical harmonics represent the deforming surfaces
* Adaptive dealiasing method is used for nonlinear operations 
* Shape reparametrization technique minimizes the high-frequency components in the spherical harmonics expansion of surface parametrization
* Weighted spherical harmonic representation treats ringing artifacts in the convection-dominated regime and provides convergent solutions
