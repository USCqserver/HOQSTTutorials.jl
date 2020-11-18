# HOQSTTutorials.jl
Tutorials for Hamiltonian Open Quantum System Toolkit(HOQST). This repo is inspired by [SciMLTutorials](https://github.com/SciML/SciMLTutorials.jl), so the workflow is the same. Because currently `HOQST` and `HOQSTTutorials.jl` are not officially registered, you need to install these packages manually.

## Table of Contents

- Introduction
  - [Introduction to HOQST through closed-system simulation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/01-closed_system.html)
  - [Introduction to HOQST through adiabatic master equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/02-single_qubit_ame.html)
  - [Redfield equation with multi-axis noise](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/03-redfield_multi_axis_noise.html)
  - [Spin-fluctuator model](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/04-spin_fluctuators.html)
  - [AME Spin-fluctuator model](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/05-ame_spin_fluctuators.html)
- Advanced
- Hamiltonian
  - [Custom eigendecomposition function](https://uscqserver.github.io/HOQSTTutorials.jl/html/hamiltonian/01-custom_eigen.html)
- Redfield
  - [Non-positivity in Redfield equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/redfield/01-non_positivity_redfield.html)
  - [Introduction to polaron transformed Redfield equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/redfield/02-polaron-transformed-redfield.html)
  - [Introduction to CGME and ULE](https://uscqserver.github.io/HOQSTTutorials.jl/html/redfield/03-CGME_ULE.html)

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