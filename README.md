# OSQATTutorials.jl
Tutorials for open system quantum annealing toolbox(OSQAT). This repo is inspired by [SciMLTutorials](https://github.com/SciML/SciMLTutorials.jl) so the workflow is the same. Because currently `OSQAT` and `OSQATTutorials.jl` are not officially registered, you need manually install these packages.

## Contributing

First of all, make sure that your current directory is `OSQATTutorials`. All
of the files are generated from the Weave.jl files in the `tutorials` folder.
To run the generation process, do for example:

```julia
using Pkg, OSQATTutorials
cd(joinpath(dirname(pathof(OSQATTutorials)), ".."))
Pkg.pkg"activate ."
Pkg.pkg"instantiate"
OSQATTutorials.weave_file("","01-closed_system.jmd")
```

To generate all of the notebooks, do:

```julia
OSQATTutorials.weave_all()
```

If you add new tutorials which require new packages, simply updating your local
environment will change the project and manifest files. When this occurs, the
updated environment files should be included in the PR.