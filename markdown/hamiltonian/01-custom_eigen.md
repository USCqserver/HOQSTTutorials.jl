---
author: "Huo Chen"
title: "Using a user-defined eigendecomposition function"
---



## Initialize a Hamiltonian with a custom eigendecomposition function

When defining a Hamiltonian object, a user-defined eigendecomposition routine can be supplied using a keyword argument `EIGS`. All the eigendecomposition calls inside the solver will call this function instead of the default one.

For example:

```julia
using OpenQuantumTools

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
There are two applications for this functionality. First, if the Hamiltonian is constant, one can precalculate that Hamiltonian's eigensystem and build a function that returns those precalculated values. This is particularly helpful for the adiabatic master equation solver.

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

Another application is to supply special eigendecomposition algorithms for sparse matrices to take advantage of the sparsity. 

For example, the default eigendecomposition algorithm for a sparse Hamiltonian is to convert it into dense matrices first and then perform the decomposition.

```julia
Hd = standard_driver(4, sp=true);
Hp = two_local_term(rand(3), [[1,2],[2,3],[3,4]], 4, sp=true)
H = SparseHamiltonian([(s)->1-s, (s)->s], [Hd, Hp], unit=:ħ)

# the default eigen_decomposition using the dense matrices algorithm
w, v = eigen_decomp(H, 0.1, lvl=4)
```

```
([-0.5740584725844259, -0.3099450687181958, -0.2975676775746419, -0.2759725
211553756], Complex{Float64}[0.23167502339168788 + 0.0im -0.000307033500147
6854 + 0.0im -0.08507953206517793 + 0.0im -0.001954564577200557 + 0.0im; -0
.24445622159282193 + 0.0im -0.19168684238925163 + 0.0im -0.1906386902185049
6 + 0.0im 0.2899748888314142 + 0.0im; … ; -0.24445622159282182 + 0.0im 0.19
168684238925138 + 0.0im 0.1906386902185051 + 0.0im -0.28997488883141437 + 0
.0im; 0.23167502339168788 + 0.0im 0.00030703350014781277 + 0.0im 0.08507953
206517775 + 0.0im 0.0019545645772008416 + 0.0im])
```




If the Hamiltonian size becomes large, we can use sparse algorithms provided by [Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) instead. Let's first load `Arpack.jl` by running:
```julia
using Arpack
```




Next, we can use an `Arpack` function to replace the default eigendecomposition routine:
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
([-0.5732983359514336, -0.3006110858203194, -0.287001096890476, -0.28629715
341853806], Complex{Float64}[-0.008572662904889262 - 0.24112608005159206im 
-0.05672192961519254 - 0.03648680658826525im -0.16320999244085255 + 0.12725
226159443936im -0.09160990155362272 - 0.06734802475623689im; 0.008861384596
505683 + 0.2492470490547713im -0.14318066811955965 - 0.09210203849372993im 
0.02376213686593718 - 0.018526964012995584im 0.34465455043365545 + 0.253376
57612672787im; … ; 0.008861384596505368 + 0.2492470490547711im 0.1431806681
195595 + 0.0921020384937301im -0.023762136865936954 + 0.018526964012995525i
m -0.34465455043365545 - 0.25337657612672765im; -0.008572662904889243 - 0.2
411260800515919im 0.056721929615192665 + 0.03648680658826529im 0.1632099924
4085263 - 0.1272522615944395im 0.09160990155362245 + 0.06734802475623684im]
)
```


