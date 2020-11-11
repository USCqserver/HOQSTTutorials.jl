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
([-0.5733365662416158, -0.29829219784009964, -0.29573131897208943, -0.27753
135589081745], Complex{Float64}[0.24006740369763377 + 0.0im -0.003646540627
8768375 + 0.0im -0.0273231196898977 + 0.0im -0.07284510111233 + 0.0im; -0.2
4909405924607872 + 0.0im -0.25786241457059506 + 0.0im -0.19831452024905785 
+ 0.0im 0.30429303692198345 + 0.0im; … ; -0.2490940592460788 + 0.0im 0.2578
6241457059483 + 0.0im 0.1983145202490578 + 0.0im -0.3042930369219833 + 0.0i
m; 0.24006740369763388 + 0.0im 0.003646540627876952 + 0.0im 0.0273231196898
97542 + 0.0im 0.07284510111233025 + 0.0im])
```




If the size of the Hamiltonian becomes very large, we can use sparse algorithms provided by [Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) instead.

```julia
using Pkg
Pkg.add("Arpack")
using Arpack

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
([-0.5734852102429372, -0.3028585932243046, -0.29341334625758725, -0.279831
9923604335], Complex{Float64}[-0.005350620073892055 - 0.23710063910047227im
 -0.004707997837175755 - 0.004206504459413104im 0.033639091854303645 + 0.09
245910660359395im -0.03174283599603078 - 0.008980832805909836im; 0.00555788
78296332936 + 0.24628524138441837im -0.14806044807728835 - 0.13228912940060
203im 0.059523708478001425 + 0.16360456255612885im 0.30263050580274736 + 0.
08562164939902227im; … ; 0.0055578878296331435 + 0.2462852413844181im 0.148
0604480772885 + 0.13228912940060195im -0.059523708478001314 - 0.16360456255
612876im -0.3026305058027473 - 0.08562164939902218im; -0.005350620073892295
 - 0.23710063910047216im 0.004707997837175793 + 0.0042065044594129415im -0.
03363909185430348 - 0.09245910660359419im 0.031742835996030755 + 0.00898083
2805909665im])
```


