<img src="assets/logo.jpg" width="256"/>

# HOQSTTutorials.jl
Tutorials for Hamiltonian Open Quantum System Toolkit(HOQST). This repo is inspired by [SciMLTutorials](https://github.com/SciML/SciMLTutorials.jl), so the workflow is the same. Because currently `HOQST` and `HOQSTTutorials.jl` are not officially registered, you need to install these packages manually.

## Table of Contents

- Introduction
  - [Introduction to HOQST through closed-system simulation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/01-closed_system.html)
  - [Introduction to HOQST through the adiabatic master equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/02-single_qubit_ame.html)
  - [Introduction to HOQST through the Lindblad equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/06-lindblad_equation.html)
  - [A tutorial on the polaron transformed Redfield equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/03-polaron_transformed_redfield.html)
  - [An Intro to the coarse-grained ME and universal Lindblad ME](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/04-CGME_ULE.html)
  - [Classical stochastic noise -- the spin-fluctuator model](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/05-spin_fluctuators.html)
- Advanced
  - Hamiltonian
    - [Custom eigendecomposition function](https://uscqserver.github.io/HOQSTTutorials.jl/html/hamiltonian/01-custom_eigen.html)
  - Redfield
    - [Non-positivity in Redfield equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/redfield/01-non_positivity_redfield.html)
    - [Redfield equation with multi-axis noise](https://uscqserver.github.io/HOQSTTutorials.jl/html/redfield/02-redfield_multi_axis_noise.html)
  - AME
    - [Adiabatic master equation with spin-fluctuators](https://uscqserver.github.io/HOQSTTutorials.jl/html/advanced/01-ame_spin_fluctuators.html)

## Contributing

First of all, make sure that your current directory is `HOQSTTutorials`. All
of the files are generated from the Weave.jl files in the `tutorials` folder.
To run the generation process, do for example:

```julia
using Pkg, HOQSTTutorials
cd(joinpath(dirname(pathof(HOQSTTutorials)), ".."))
Pkg.pkg"activate ."
Pkg.pkg"instantiate"
HOQSTTutorials.weave_file("introduction","01-closed_system.jmd")
```

To generate all of the notebooks, do:

```julia
HOQSTTutorials.weave_all()
```

If you add new tutorials that require new packages, simply updating your local
environment will change the project and manifest files. When this occurs, the
updated environment files should be included in the PR.