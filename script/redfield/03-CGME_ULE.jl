
using OrdinaryDiffEq, Plots, LaTeXStrings
using QuantumAnnealingTools

# Hamiltonian
H = DenseHamiltonian([(s)->1-s, (s)->s], -[σx, σz]/2, unit=:ħ)
# initial state
u0 = PauliVec[1][1]
# coupling
coupling = ConstantCouplings(["Z"], unit=:ħ)
# bath
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)

tf = 60
U = solve_unitary(annealing, tf, alg=Tsit5(), abstol=1e-8, reltol=1e-8)
U = InplaceUnitary(U)

@time solr = solve_redfield(annealing, tf, U, alg=Tsit5())
# we set the integration error tolerance to 1e-5 for speed
@time solc = solve_CGME(annealing, tf, U, alg=Tsit5(), int_atol=1e-5, int_rtol=1e-5)
plot(solr, H, [1], 0:0.01:tf, linewidth=2, xlabel="t (ns)", ylabel="\$P_G(t)\$", label="Redfield")
plot!(solc, H, [1], 0:0.01:tf, linewidth=2, label="CGME")


using QuadGK

g(t) = quadgk((w)->sqrt(γ(w, bath))*exp(1.0im*w*t)/2/π, -Inf, Inf)[1]

t = range(-0.5,0.5,length=500)
g_value = g.(t)
c_value = [correlation(x, bath) for x in t];

plot(t, real.(g_value), label="Re[g(t)]", linewidth=2)
plot!(t, imag.(g_value), label="Im[g(t)]", linewidth=2)
plot!(t, real.(c_value), label="Re[C(t)]", linewidth=2)
plot!(t, imag.(c_value), label="Im[C(t)]", linewidth=2)
xlabel!("t (ns)")
ylabel!("correlation")


t = range(-4,4,length=2000)
g_value = g.(t)
gf = construct_interpolations(t, g_value, extrapolation = "flat")

t = range(-0.5,0.5,length=500)
g_value = gf.(t)
plot(t, real.(g_value), label="Re[g(t)]", linewidth=2)
plot!(t, imag.(g_value), label="Im[g(t)]", linewidth=2)
xlabel!("t (ns)")
ylabel!("fitted jump correlator")


ubath = ULEBath(gf)
annealing = Annealing(H, u0; coupling=coupling, bath=ubath)
@time solu = solve_ule(annealing, tf, U, alg=Tsit5(), int_atol=1e-5, int_rtol=1e-5)
plot(solr, H, [1], 0:0.01:tf, linewidth=2, xlabel="t (ns)", ylabel="\$P_G(t)\$", label="Redfield")
plot!(solc, H, [1], 0:0.01:tf, linewidth=2, label="CGME")
plot!(solu, H, [1], 0:0.01:tf, linewidth=2, label="ULE")

