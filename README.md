<div style="padding:2px;background:#fff;webkit-border-radius:10px;-moz-border-radius:10px;border-radius:10px;width:60%;margin:0 auto;overflow:hidden;">
<p align="center">
  <a href="http://mflowcode.github.io/">
    <img src="install/cells-2.png" alt="RBC3D Banner" width="500"/>
  </a>
</p>
</div>

# Spectral Boundary Integral Solver for Cell-scale Flows
## Authors: S. H. Bryngelson, H. Zhao, A. Isfahani, J. B. Freund

[![DOI](https://zenodo.org/badge/412637841.svg)](https://zenodo.org/badge/latestdoi/412637841)


RBC3D is a flow solver for soft capsules and cells. 
It implements the methods discussed in [Zhao et al., JCP (2010)](https://doi.org/10.1016/j.jcp.2010.01.024) and more.
In particular, it solves the boundary integral form of the Stokes equations via an algorithm tailored for cell-scale simulations:

* Spectrally-accurate spherical harmonics represent the deforming surfaces
* Modified Green’s function approximation used for near-range interactions
* Electrostatic-like repulsion prevents cells from intersecting
* Weak-formulation of no-slip boundary conditions (e.g., vessel walls)
* These features ensure that simulations are robust. Parallel communication (MPI) enables large simulations, such as model vascular networks.
