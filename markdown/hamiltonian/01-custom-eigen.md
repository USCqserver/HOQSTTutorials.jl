---
author: "Huo Chen"
title: "Using a user defined eigendecomposition function"
---


## Initialize Hamiltonian with a custom eigendecomposition function
When defining a Hamiltonian object, a keyword argument `EIGS` can be used to construct a user defined eigendecomposition routine. All the eigendecomposition calls inside the solver will call this function instead of the default one.

For example:

````julia
using QuantumAnnealingTools

# define a function to construct an eigendecomposition routine for the
# Hamiltonian. The routine should have signature: (H, t, lvl) -> (w, v).
# the argument is the cache of the Hamiltonian
function build_user_eigen(u_cache)
    EIGS = function(H, t, lvl)
        println("I am the user defined eigendecomposition routine.")
        w, v = eigen(Hermitian(H(t)))
        w[1:lvl], v[:, 1:lvl]
    end
end

H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz], EIGS=build_user_eigen)

eigen_decomp(H, 0.5, lvl=2)
````


````
I am the user defined eigendecomposition routine.
([-0.7071067811865476, 0.7071067811865476], Complex{Float64}[-0.38268343236
508984 + 0.0im 0.9238795325112866 - 0.0im; 0.9238795325112868 - 0.0im 0.382
6834323650897 - 0.0im])
````





### Constant Hamiltonian
There are two applications for this functionality. First, if the Hamiltonian is constant, one can precalculates the eigensystem of that Hamiltonian and build a function that returns those precalculated values. This is particular helpful for adiabatic master equation solver.

````julia
function build_user_eigen(u_cache)
    # note that to keep the unit consistent, the unit of any value inside the routine should be 1/h
    w, v = eigen(Hermitian(2*π*(σx+σz)))
    EIGS = function(H, t, lvl)
        w[1:lvl], v[:, 1:lvl]
    end
end

H = DenseHamiltonian([(s)->1.0], [σx+σz], EIGS=build_user_eigen)

print(eigen_decomp(H, 0.5, lvl=2))
````


````
([-1.4142135623730951, 1.4142135623730951], Complex{Float64}[0.382683432365
0897 + 0.0im -0.9238795325112867 + 0.0im; -0.9238795325112867 + 0.0im -0.38
26834323650897 + 0.0im])
````



````julia
print(eigen_decomp(H, 0.0, lvl=2))
````


````
([-1.4142135623730951, 1.4142135623730951], Complex{Float64}[0.382683432365
0897 + 0.0im -0.9238795325112867 + 0.0im; -0.9238795325112867 + 0.0im -0.38
26834323650897 + 0.0im])
````





### Sparse Hamiltonian

Another application is to supply speical eigendecomposition algorithms for sparse matrices in order to take advantage of the sparsity. 

For example, the default eigendecomposition algorithm for sparse Hamiltonian is to first convert it into dense matrices and then perform the decomposition.

````julia
Hd = standard_driver(4, sp=true);
Hp = two_local_term(rand(3), [[1,2],[2,3],[3,4]], 4, sp=true)
H = SparseHamiltonian([(s)->1-s, (s)->s], [Hd, Hp], unit=:ħ)

# the default eigen_decomposition using dense matrices algorithm
w, v = eigen_decomp(H, 0.1, lvl=4)
````


````
([-0.5737621770602794, -0.3069723449558068, -0.29469862289106913, -0.278699
86805107744], Complex{Float64}[0.2341949813011429 + 0.0im 0.004930130285252
878 + 0.0im -0.10402072727232821 + 0.0im -0.024394560895135312 + 0.0im; -0.
24404906262614626 + 0.0im -0.1728822289031417 + 0.0im -0.18971279526222828 
+ 0.0im -0.27744707653185785 + 0.0im; … ; -0.24404906262614612 + 0.0im 0.17
288222890314184 + 0.0im 0.18971279526222817 + 0.0im 0.2774470765318579 + 0.
0im; 0.23419498130114297 + 0.0im -0.004930130285252922 + 0.0im 0.1040207272
7232832 + 0.0im 0.024394560895135323 + 0.0im])
````




If the size of the Hamiltonian becomes very large, we can use sparse algorithms provided by [Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) instead.

````julia
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
````


````
Using sparse matrix algorithm
([-0.5733249332272754, -0.3008323040413639, -0.2898582881644965, -0.2833998
0361983713], Complex{Float64}[-0.23712984825025823 - 0.03491516672453986im 
-0.02506183488516595 + 0.02496150965046808im -0.12985402694859677 - 0.03566
199745553825im -0.10176624397424387 + 0.06486367296562603im; 0.240860443966
6072 + 0.03546446227876591im 0.08277783806062194 - 0.08244646942501556im -0
.19576059316961733 - 0.05376201215750451im -0.17491925665628236 + 0.1114898
7145497774im; … ; 0.2408604439666074 + 0.03546446227876591im -0.08277783806
062196 + 0.08244646942501518im 0.19576059316961733 + 0.053762012157504535im
 0.1749192566562824 - 0.11148987145497787im; -0.2371298482502583 - 0.034915
166724539905im 0.02506183488516586 - 0.02496150965046839im 0.12985402694859
69 + 0.03566199745553817im 0.10176624397424397 - 0.06486367296562622im])
````


