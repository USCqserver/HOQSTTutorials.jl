
using OrdinaryDiffEq, QuantumAnnealingTools, Plots
using LaTeXStrings

function err_bound(tf, cfun)
    tsb, esb = τ_SB(cfun)
    tb, eb = τ_B(cfun, tf, tsb)
    tb / tsb
end

fc = 4; T =12; tf = 1000;
ηlist = log_uniform(1e-3, 5, 1000)
err_ratio = []
err_clist = []
err_klist = []
for η in ηlist
    bath = Ohmic(η, fc, T)
    cfun = (x)->correlation(x, bath)
    pfun = (x)->polaron_correlation(x, bath)
    err_c = err_bound(tf, cfun)
    err_k = err_bound(tf, pfun)
    push!(err_clist, err_c)
    push!(err_klist, err_k)
    push!(err_ratio, err_c/err_k)
end
idx = findfirst((x)->x>=1, err_ratio)
plot(ηlist, err_ratio, xscale=:log10, yscale=:log10, label="", linewidth=2)
vline!([ηlist[idx]], label="", linestyle=:dash, linewidth=2)
annotate!([(0.5, 1.0, Plots.text("polaron")), (0.01, 1.0, Plots.text("Redfield"))])
xlabel!(L"\eta g^2")
ylabel!("R")


plot(ηlist, err_clist, xscale=:log10, yscale=:log10, label="Redfield", linewidth=2)
plot!(ηlist, err_klist, xscale=:log10, yscale=:log10, label="PTRE", linewidth=2)
xlabel!(L"\eta g^2")
ylabel!("error")


    # assume ϵ = 1
    const Δ = 0.1 
    # define the Ohmic bath in polaron transformed frame
    η = 0.5; bath = Ohmic(η, fc, T)
    K(t1, t2) = Δ^2 * polaron_correlation(t1-t2, bath)
    cfun = [nothing K; K nothing]
    pbath = CorrelatedBath(((1,2),(2,1)), correlation=cfun)
    # define coupling as σ+ and σ- operators
    σp = [0 1;0 0.0im]; σm = [0 0;1 0.0im]
    coupling = ConstantCouplings([σp, σm])
    # manually define the unitary operator
    U(t) = exp(-2.0im * π * σz * t)
    H = DenseHamiltonian([(s)->1.0], [σz])
    u0 = PauliVec[3][1]
    annealing = Annealing(H, u0, coupling = coupling, bath = pbath)
    tf = 100
    sol_ptre = solve_redfield(annealing, tf, U, alg=Tsit5(), Ta=2, reltol=1e-5)
    pop_e = [real(s[1,1]) for s in sol_ptre.u]
    plot(sol_ptre.t, pop_e, xlabel=L"t\ (\mathrm{ns})", ylabel=L"P_0(t)", label="", linewidth = 2)


H = DenseHamiltonian([(s)->1.0], [σz+0.1*σx])
coupling = ConstantCouplings(["Z"])
annealing = Annealing(H, u0, coupling = coupling, bath = bath)
tf = 100
sol_redfield = solve_redfield(annealing, tf, U, alg=Tsit5(), Ta=40, reltol=1e-5, callback=PositivityCheckCallback())
pop_e_redfield = [real(s[1,1]) for s in sol_redfield.u]
plot(sol_ptre.t, pop_e, xlabel=L"t\ (\mathrm{ns})", ylabel=L"P_0(t)", label="PTRE", linewidth = 2, legend = :right)
plot!(sol_redfield.t, pop_e_redfield, xlabel=L"t\ (\mathrm{ns})", ylabel=L"P_0(t)", label="Redfield", linewidth = 2)


t_axis = range(0, 5, length=100)
off_diag_ptre = [abs(sol_ptre(t)[1,2]) for t in t_axis]
off_diag_redfield = [abs(sol_redfield(t)[1,2]) for t in t_axis]
plot(t_axis, off_diag_ptre, xlabel=L"t\ (\mathrm{ns})", ylabel=L"\lvert \rho_{12} \rvert|(t)", label="PTRE", linewidth = 2, legend=:right)
plot!(t_axis, off_diag_redfield, xlabel=L"t\ (\mathrm{ns})", ylabel=L"|\rho_{12}(t)|", label="Redfield", linewidth = 2)

