
using OrdinaryDiffEq, QuantumAnnealingTools, Plots

coupling_1 = ConstantCouplings(["X"])
bath_1 = Ohmic(1e-4, 4, 16)
interaction_1 = Interaction(coupling_1, bath_1)

coupling_2 = ConstantCouplings(["Z"])
bath_2 = Ohmic(1e-4, 0.1, 16)
interaction_2 = Interaction(coupling_2, bath_2)

interaction_set = InteractionSet(interaction_1, interaction_2);


H = DenseHamiltonian([(s) -> 1.0], -[σz], unit = :ħ)
u0 = PauliVec[1][1]
annealing_1 = Annealing(H, u0, coupling=coupling_1, bath = bath_1)
annealing_2 = Annealing(H, u0, coupling=coupling_2, bath = bath_2)
annealing = Annealing(H, u0, interactions=interaction_set)


tf = 10
# Generate the unitary first
U = solve_unitary(annealing, tf, alg=Tsit5(), reltol=1e-6)
# tag the unitary so the solver know it has inplace update method
# this will speed up the calculation of integral
U = InplaceUnitary(U)
 # Solve the Redfield equation
sol_1 = solve_redfield(annealing_1, tf, U, alg = Tsit5(), abstol = 1e-6, reltol = 1e-6)
sol_2 = solve_redfield(annealing_2, tf, U, alg = Tsit5(), abstol = 1e-6, reltol = 1e-6)
sol = solve_redfield(annealing, tf, U, alg = Tsit5(), abstol = 1e-6, reltol = 1e-6)


t_list = range(0,tf,length=200)
x1 = []
x2 = []
x = []
for s in t_list
    push!(x1, real(tr(σx*sol_1(s))))
    push!(x2, real(tr(σx*sol_2(s))))
    push!(x, real(tr(σx*sol(s))))
end
x_nopulse = x1
plot(t_list, x1, linewidth=2, label="X coupling")
plot!(t_list, x2, linewidth=2, label="Z coupling")
plot!(t_list, x, linewidth=2, label="X+Z coupling")
xlabel!("s")
ylabel!("<X>")


# in this example, we apply an Z pulse in the middle of the annealing
# for the InstPulseCallback constructor
# the first argument is a list of times where the pulses are applied
# the second argument is a function to update the state update!(c, pulse_index
# the function will update the state c with give pulse_index
cbu = InstPulseCallback([0.5 * tf], (c, x) -> c .= σx * c)
cb = InstPulseCallback([0.5 * tf], (c, x) -> c .= σx * c * σx)
annealing_1 = Annealing(H, u0, coupling = coupling_1, bath = bath_1)
annealing_2 = Annealing(H, u0, coupling=coupling_2, bath = bath_2)
annealing = Annealing(H, u0, interactions=interaction_set)

tf = 10
U = solve_unitary(annealing, tf, alg=Tsit5(), reltol=1e-6, callback = cbu);
U = InplaceUnitary(U)


sol_1 = solve_redfield(annealing_1, tf, U, alg = Tsit5(), reltol = 1e-6, callback=cb)
sol_2 = solve_redfield(annealing_2, tf, U, alg = Tsit5(), reltol = 1e-6, callback=cb)
sol = solve_redfield(annealing, tf, U, alg = Tsit5(), reltol = 1e-6, callback=cb);
t_list = range(0,tf,length=200)
x1 = []
x2 = []
x = []
for s in t_list
    push!(x1, real(tr(σx*sol_1(s))))
    push!(x2, real(tr(σx*sol_2(s))))
    push!(x, real(tr(σx*sol(s))))
end
plot(t_list, x1, linewidth=2, label="X coupling", legend=:topleft)
plot!(t_list, x2, linewidth=2, label="Z coupling")
plot!(t_list, x, linewidth=2, label="X+Z coupling")
plot!(t_list, x_nopulse, linewidth=2, linestyle=:dash, label="X coupling with no pulse")
xlabel!("s")
ylabel!("<X>")

