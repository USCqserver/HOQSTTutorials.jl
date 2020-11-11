---
author: "Huo Chen"
title: "Using a user-defined eigendecomposition function"
---


## Initialize Hamiltonian with a custom eigendecomposition function
When defining a Hamiltonian object, a keyword argument `EIGS` can be supplied to construct a user defined eigendecomposition routine. All the eigendecomposition calls inside the solver will call this function instead of the default one.

For example:

```julia
using QuantumAnnealingTools

# Define a function to construct an eigendecomposition routine for the
# Hamiltonian. The routine should have the signature: (H, t, lvl) -> (w, v).
# The argument of this function is the cache used by the Hamiltonian object.
function build_user_eigen(u_cache)
    EIGS = function(H, t, lvl)
        println("I am the user defined eigendecomposition routine.")
        w, v = eigen(Hermitian(H(t)))
        w[1:lvl], v[:, 1:lvl]
    end
end

H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz], EIGS=build_user_eigen)

eigen_decomp(H, 0.5, lvl=2)
```

```
I am the user defined eigendecomposition routine.
([-0.7071067811865476, 0.7071067811865476], Complex{Float64}[-0.38268343236
508984 + 0.0im 0.9238795325112866 - 0.0im; 0.9238795325112868 - 0.0im 0.382
6834323650897 - 0.0im])
```





### Constant Hamiltonian
There are two applications for this functionality. First, if the Hamiltonian is constant, one can precalculates the eigensystem of that Hamiltonian and build a function that returns those precalculated values. This is particular helpful for adiabatic master equation solver.

```julia

function build_user_eigen(u_cache)
    # note that to keep the unit consistent, the unit of any value inside the routine should be 1/h
    w, v = eigen(Hermitian(2*π*(σx+σz)))
    EIGS = function(H, t, lvl)
        w[1:lvl], v[:, 1:lvl]
    end
end

H = DenseHamiltonian([(s)->1.0], [σx+σz], EIGS=build_user_eigen)

print(eigen_decomp(H, 0.5, lvl=2))
print(eigen_decomp(H, 0.0, lvl=2))
```

```
([-1.4142135623730951, 1.4142135623730951], Complex{Float64}[0.382683432365
0897 + 0.0im -0.9238795325112867 + 0.0im; -0.9238795325112867 + 0.0im -0.38
26834323650897 + 0.0im])([-1.4142135623730951, 1.4142135623730951], Complex
{Float64}[0.3826834323650897 + 0.0im -0.9238795325112867 + 0.0im; -0.923879
5325112867 + 0.0im -0.3826834323650897 + 0.0im])
```





### Sparse Hamiltonian

Another application is to supply speical eigendecomposition algorithms for sparse matrices in order to take advantage of the sparsity. 

For example, the default eigendecomposition algorithm for sparse Hamiltonian is to first convert it into dense matrices and then perform the decomposition.

```julia
Hd = standard_driver(4, sp=true);
Hp = two_local_term(rand(3), [[1,2],[2,3],[3,4]], 4, sp=true)
H = SparseHamiltonian([(s)->1-s, (s)->s], [Hd, Hp], unit=:ħ)

# the default eigen_decomposition using dense matrices algorithm
w, v = eigen_decomp(H, 0.1, lvl=4)
```

```
([-0.5729923320369654, -0.2909140176773721, -0.28692493464359575, -0.286089
30359797324], Complex{Float64}[0.2473784274380235 + 0.0im 0.040695944338657
37 + 0.0im 0.05915731187328016 + 0.0im -0.28180804014959693 + 0.0im; -0.251
0321427724954 + 0.0im 0.2932926473618397 + 0.0im 0.04998787555385241 + 0.0i
m 0.3971501566526018 + 0.0im; … ; -0.2510321427724953 + 0.0im -0.2932926473
618393 + 0.0im -0.04998787555385233 + 0.0im -0.3971501566526014 + 0.0im; 0.
24737842743802355 + 0.0im -0.04069594433865745 + 0.0im -0.05915731187327999
4 + 0.0im 0.28180804014959704 + 0.0im])
```




If the size of the Hamiltonian becomes very large, we can use sparse algorithms provided by [Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) instead. Let's first install `Arpack.jl` by running:
```julia
using Pkg
Pkg.add("Arpack")
using Arpack
```




Next, we can use `Arpack` function to replace the default eigendecomposition routine:
```julia
function build_user_eigen(u_cache)
    function (H, t, lvl)
        hmat = H(t)
        println("Using sparse matrix algorithm")
        # define all the Arpack routine parameters here
        eigs(hmat, nev = lvl, which=:SR, tol=0.0, maxiter=300)
    end
end

Hd = standard_driver(4, sp=true);
Hp = two_local_term(rand(3), [[1,2],[2,3],[3,4]], 4, sp=true)
H = SparseHamiltonian([(s)->1-s, (s)->s], [Hd, Hp], unit=:ħ, EIGS =build_user_eigen)

eigen_decomp(H, 0.1, lvl=4)
```

```
Using sparse matrix algorithm
([-0.5733165321046936, -0.3004703651326788, -0.29063522075372417, -0.282563
26509951385], Complex{Float64}[-0.007159439655548368 - 0.23936480997673287i
m -0.012133557946698244 + 0.007450789819478774im 0.02035383149534172 - 0.13
33297763551637im -0.034271552359772174 - 0.03703217538616028im; 0.007373209
090542359 + 0.24651186095390382im -0.15570674055672995 + 0.0956140154817683
2im 0.021833559400207706 - 0.14302287962511612im 0.23063593297888502 + 0.24
92139903896338im; … ; 0.007373209090542037 + 0.24651186095390368im 0.155706
74055673003 - 0.09561401548176823im -0.02183355940020778 + 0.14302287962511
614im -0.2306359329788852 - 0.2492139903896334im; -0.00715943965554853 - 0.
23936480997673287im 0.0121335579466981 - 0.007450789819478828im -0.02035383
1495342048 + 0.1333297763551638im 0.03427155235977204 + 0.03703217538615996
6im])
```


