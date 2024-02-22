# Trixi simulations

This directory contains code and instructions to reproduce the numerical
experiments reported in the article

> An Entropy-Stable Discontinuous Galerkin Discretization of the Ideal Multi-Ion Magnetohydrodynamics System

The results were obtained using Julia v1.10.0 and [this version of Trixi.jl](https://github.com/amrueda/Trixi.jl/tree/9e1dbfff4fcd3b80b2692c39950efa93d76e8a04). When this reproducibility repository was last updated, the implementation of the multi-ion GLM-MHD equations wa s not yet merged into the main [Trixi.jl](https://github.com/trixi-framework/Trixi.jl/) repository.

## Instructions

To run the examples, follow the instructions:

* Move to this directory and clone the fork of Trixi.jl repository:
  ```bash
  git clone git@github.com:amrueda/Trixi.jl.git
  ```
* Move to the Trixi.jl folder and change to the branch with which the examples were computed:
  ```bash
  cd Trixi.jl
  git checkout 9e1dbfff4fcd3b80b2692c39950efa93d76e8a04
  ```
* Create a run direcrtoty and install all the dependencies of Trixi.jl:
  ```bash
  mkdir run
  cd run
  julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path=".."))' # Install local Trixi.jl clone
  julia --project=. -e 'using Pkg; Pkg.add(["OrdinaryDiffEq", "Trixi2Vtk", "Plots"])' # Install additional packages
  ```
* Run the individual examples using julia, e.g.,
  ```bash
  julia --check-bounds=no --threads=10 -e 'include(joinpath("..", "convergence", "runelixirs.jl"))'
  ```
* The simulation files will be stored in separate output folders. To visualize with paraview use Trixi2Vtk, e.g.,:
  ```bash
  julia --check-bounds=no --threads=10 -e 'using Trixi2Vtk ; trixi2vtk(joinpath("..", "khi", "khi_glm_H_H2", "out", "solution_0*.h5"), output_directory="out")'
  ```
