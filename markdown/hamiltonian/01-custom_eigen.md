---
author: "Huo Chen"
title: "Using a user-defined eigendecomposition function"
---



## Initialize Hamiltonian with a custom eigendecomposition function

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

For example, the default eigendecomposition algorithm for sparse Hamiltonian is to convert it into dense matrices first and then perform the decomposition.

```julia
Hd = standard_driver(4, sp=true);
Hp = two_local_term(rand(3), [[1,2],[2,3],[3,4]], 4, sp=true)
H = SparseHamiltonian([(s)->1-s, (s)->s], [Hd, Hp], unit=:ħ)

# the default eigen_decomposition using dense matrices algorithm
w, v = eigen_decomp(H, 0.1, lvl=4)
```

```
([-0.5735244654671962, -0.304565310552158, -0.28936834479657075, -0.2842127
439114835], Complex{Float64}[0.23764845654771619 + 0.0im -0.009503985289811
945 + 0.0im 0.19830147707586077 + 0.0im -0.019056568660583945 + 0.0im; -0.2
4261757772878467 + 0.0im 0.11756923861584401 + 0.0im 0.12692646580471256 + 
0.0im -0.30325080327059184 + 0.0im; … ; -0.2426175777287845 + 0.0im -0.1175
6923861584388 + 0.0im -0.12692646580471267 + 0.0im 0.3032508032705917 + 0.0
im; 0.23764845654771624 + 0.0im 0.00950398528981197 + 0.0im -0.198301477075
86065 + 0.0im 0.019056568660583927 + 0.0im])
```




If the Hamiltonian size becomes large, we can use sparse algorithms provided by [Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) instead. Let's first install `Arpack.jl` by running:
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
([-0.5736067127748863, -0.305375099581338, -0.29218673198844436, -0.2812376
250417108], Complex{Float64}[0.0077968247914967405 - 0.2359246876933563im 0
.01778350602287372 + 0.018094339548598005im -0.12077612792231388 + 0.011769
929520151275im -0.01276513438898149 - 0.09912659311490553im; -0.00800540859
5379472 + 0.24223624016566098im -0.09281537302257294 - 0.09443767008824275i
m -0.20280559312074686 + 0.019763901843741827im -0.02829699560284785 - 0.21
973797408028825im; … ; -0.008005408595379383 + 0.24223624016566106im 0.0928
1537302257299 + 0.09443767008824264im 0.20280559312074678 - 0.0197639018437
41636im 0.028296995602847912 + 0.21973797408028825im; 0.007796824791496877 
- 0.23592468769335637im -0.017783506022873815 - 0.01809433954859785im 0.120
77612792231385 - 0.0117699295201512im 0.012765134388981386 + 0.099126593114
90575im])
```


