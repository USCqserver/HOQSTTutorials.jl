
using OpenQuantumTools, OrdinaryDiffEq, Plots
# define the Hamiltonian
H = DenseHamiltonian([(s)->1.0], [σz], unit=:ħ)
# define the initial state
u0 = PauliVec[1][1]*PauliVec[1][1]'
# define the Lindblad operator
# the rate and Lindblad operator can also be time-dependent functions
lind = Lindblad(0.1, σz)
# combine them into an Annealing object
annealing = Annealing(H, u0, interactions = InteractionSet(lind))


# define total annealing/evolution time
tf = 10
# solve the Lindblad equation
sol = solve_lindblad(annealing, 10, alg=Tsit5());


t_axis = range(0, 10, length=100)
bloch_vector = []
for t in t_axis
    # matrix_decompose projects a matrix onto a list of basis elements
    push!(bloch_vector, 2*real.(matrix_decompose(sol(t), [σx, σy, σz])))
end

off_diag = []
for t in t_axis
    push!(off_diag, abs(sol(t)[1,2]))
end


plot(t_axis, [c[1] for c in bloch_vector], label="X", linewidth=2)
plot!(t_axis, [c[2] for c in bloch_vector], label="Y", linewidth=2)
plot!(t_axis, [c[3] for c in bloch_vector], label="Z", linewidth=2)
xlabel!("t (ns)")
ylabel!("Bloch Vector")


plot(t_axis, off_diag, linewidth=2, label="ME")
plot!(t_axis, 0.5*exp.(-0.2*t_axis), linestyle=:dash, linewidth=3, label="Analytical")
xlabel!("t (ns)")
ylabel!("|ρ₀₁(t)|")


# For the quantum trajectories method, the u0 supplied to `Annealing` must be
# a state vector.
# We will show how to replace it with a pure state ensemble later
# in this example
u0 = PauliVec[1][1]
lind = Lindblad(0.1, σz)
annealing = Annealing(H, u0, interactions = InteractionSet(lind))
tf = 10
prob = build_ensembles(annealing, tf, :lindblad)
# We ran each trajectory serially for the sake of simplicity. The user is encouraged to try parallel algorithms.
sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1000, saveat=range(0,tf,length=100))


vec_norm = []
# this is the index of the trajectory you want to look at
idx = 2
for v in sol[idx].u
    push!(vec_norm, norm(v))
end
plot(sol[idx].t, vec_norm, linewidth=2, label="", xlabel="t (ns)", ylabel="‖ψ̃(t)‖")


t_axis = range(0,tf,length=100)
dataset = []
for t in t_axis
    temp = []
    for so in sol
        v = normalize(so(t))
        push!(temp, real.(v'*σx*v))
    end
    push!(dataset, temp)
end

x_mean = []
x_sem = []
for data in dataset
    t_mean = sum(data)/1000
    t_sem = sqrt(sum((x)->(x-t_mean)^2, data))/1000
    push!(x_mean, t_mean)
    push!(x_sem, t_sem)
end

scatter(t_axis, x_mean, marker=:d, yerror=2 * x_sem, label="1000 trajectories", markersize=6, ylabel="<X>", xlabel="t (ns)")
plot!(t_axis, [c[1] for c in bloch_vector], linewidth=2, label="direct solver")

dataset = []
for t in t_axis
    temp = []
    for so in sol[1:100]
        v = normalize(so(t))
        push!(temp, real.(v'*σx*v))
    end
    push!(dataset, temp)
end

x_mean = []
x_sem = []
for data in dataset
    t_mean = sum(data)/100
    t_sem = sqrt(sum((x)->(x-t_mean)^2, data))/100
    push!(x_mean, t_mean)
    push!(x_sem, t_sem)
end

plot!(t_axis, x_mean, ribbon=2*x_sem, label="100 trajectories", markersize=6, ylabel="<X>", xlabel="t (ns)")


# PuliVec[1][1] is the plus state and PauliVec[1][2] is the minus state
# The first argument is a list of corresponding probabilities of the
# pure states in the second argument. 
E = EᵨEnsemble([0.7, 0.3], [PauliVec[1][1], PauliVec[1][2]])

function prob_func(prob,i,repeat)
  prob.u0 .= sample_state_vector(E)
  prob
end

u0 = PauliVec[1][1]
lind = Lindblad(0.1, σz)
annealing = Annealing(H, u0, interactions = InteractionSet(lind))
tf = 10
prob = build_ensembles(annealing, tf, :lindblad, prob_func=prob_func)
sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=2000, saveat=range(0,tf,length=100))


initial_state_counter = zeros(2)
for so in sol
    if so.prob.u0 == PauliVec[1][1]
        initial_state_counter[1] += 1
    else
        initial_state_counter[2] += 1
    end
end
bar([0,1],initial_state_counter, label="", ylabel="Frequency", xticks=([0, 1], ["|+⟩","|-⟩"]))


t_axis = range(0,tf,length=100)
dataset = []
for t in t_axis
    temp = []
    for so in sol
        v = normalize(so(t))
        push!(temp, real.(v'*σx*v))
    end
    push!(dataset, temp)
end

x_mean = []
x_sem = []
for data in dataset
    t_mean = sum(data)/2000
    t_sem = sqrt(sum((x)->(x-t_mean)^2, data))/2000
    push!(x_mean, t_mean)
    push!(x_sem, t_sem)
end

# define Hamiltoian
H = DenseHamiltonian([(s)->1.0], [σz], unit=:ħ)
# define initial state
u0 = 0.7*PauliVec[1][1]*PauliVec[1][1]'+0.3*PauliVec[1][2]*PauliVec[1][2]'
# define Lindblad operator
lind = Lindblad(0.1, σz)
# combine them into an Annealing object
annealing = Annealing(H, u0, interactions = InteractionSet(lind))

# define total annealing/evolution time
tf = 10
# solve the Lindblad equation
sol = solve_lindblad(annealing, 10, alg=Tsit5());

t_axis = range(0, 10, length=100)
x_vector = []
for t in t_axis
    push!(x_vector, 2*real.(matrix_decompose(sol(t), [σx])))
end

scatter(t_axis, x_mean, marker=:d, yerror=2*x_sem, label="2000 trajectories", markersize=6, ylabel="<X>", xlabel="t (ns)")

plot!(t_axis, [c[1] for c in x_vector], linewidth=2, label="direct solver")


using HOQSTTutorials
HOQSTTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

