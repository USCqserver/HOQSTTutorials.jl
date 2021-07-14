
using OpenQuantumTools, LaTeXStrings
using OrdinaryDiffEq, Plots, StatsBase

# calculate the mean and standard deviation of the mean estimator from a sample
function mean_std(sample)
    m, v = mean_and_std(sample, corrected=true)
    m, v/sqrt(length(sample))
end

# All values calculated below are in angular frequency units
num = 10
bvec = 0.2 * ones(num)
γvec = log_uniform(0.01, 1, num)
fluctuator_ensemble = EnsembleFluctuator(bvec, γvec);

plot(fluctuator_ensemble, :spectrum, 2*π*range(0.01, 4, length=100), xscale=:log10, yscale=:log10, linewidth=2, label="")
xlabel!(L"\omega")
ylabel!(L"S(\omega)")


H = DenseHamiltonian([(s)->1], [σz], unit=:ħ)
coupling = ConstantCouplings([0.5*σz], unit=:ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0, coupling = coupling, bath=fluctuator_ensemble)
tf = 10
# create object for parallel simulation
prob = build_ensembles(annealing, tf, :stochastic)
# we run each trajectory serially for this example
sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1000, reltol=1e-6, saveat=range(0,tf,length=100))


t_axis = range(0,tf,length=100)
es = []
err = []
for (i, s) in enumerate(t_axis)
    sample = [abs2(u0'*so[i]) for so in sol]
    pop, pop_std = mean_std(sample)
    push!(es, pop)
    push!(err, 2*pop_std)
end

plot(t_axis, es, ribbon=err, linewidth=2, label="")
xlabel!("ns")
ylabel!("<+>")


cb = InstPulseCallback([0.25, 0.5, 0.75] .* tf, (c, x) -> c .= σx * c)
sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1000, reltol=1e-6, saveat=range(0,tf,length=100), callback=cb)


s_axis = range(0,1,length=100)
es = []
err = []
for (i, s) in enumerate(s_axis)
    sample = [abs2(u0'*so[i]) for so in sol]
    pop, pop_std = mean_std(sample)
    push!(es, pop)
    push!(err, 2*pop_std)
end

plot(tf*s_axis, es, ribbon=err, linewidth=2, label="")
xlabel!("ns")
ylabel!("<+>")


using HOQSTTutorials
HOQSTTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

