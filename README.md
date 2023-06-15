# Spectral Boundary Integral Solver for Electrohydrodynamic Flows in Viscous Drops
## Authors: M. Firouznia, S. H. Bryngelson, D. Saintillan

EHD_Drop_3D is an electrohydrodynamic solver for viscous drops. It solves the boundary integral forms of the Laplace's equations and Stokes equations via an algorithm tailored for drop-scale simulations:

* Spectrally-accurate spherical harmonics represent the deforming surfaces
* Adaptive dealiasing method is used for nonlinear operations 
* Shape reparametrization technique minimizes the high-frequency components in the spherical harmonics expansion of surface parametrization
* Weighted spherical harmonic representation treats ringing artifacts in the convection-dominated regime and provides convergent solutions

For more details about EHD_Drop_3D please refer to the associated paper in the Journal of Computational Physics ([Firouznia et al., JCP (2023)](https://doi.org/10.1016/j.jcp.2023.112248)). This code builds off of [RBC3D](https://github.com/comp-physics/RBC3D), extending it to model interfacial charge transport in viscous drops. RBC3D was written by H. Zhao, A. Isfahani, S. H. Bryngelson, and J. B. Freund, which has an associated JCP paper: [Zhao et al., JCP (2010)](https://doi.org/10.1016/j.jcp.2010.01.024).
