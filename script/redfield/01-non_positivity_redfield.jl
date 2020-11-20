
using OpenQuantumTools, OrdinaryDiffEq, Plots, Printf, LaTeXStrings

β = 4
T = β_2_temperature(β)
η = 0.1
fc= 10/(2π)
bath = Ohmic(η, fc, T)


plot(bath, :γ, range(0,10,length=100), linewidth=2, label="")


τsb, err_τsb = τ_SB((x)->correlation(x, bath))
@printf("τ_sb of the Ohmic bath is %.6f with error estimation %.2e \n", τsb, err_τsb)
τb, err_τb = τ_B((x)->correlation(x, bath), 100, τsb)
@printf("τ_b of the Ohmic bath is %.6f with error estimation %.2e \n", τb, err_τb)


Hp = 0.5*σz⊗σi - 0.7*σi⊗σz + 0.3*σz⊗σz
Hd = standard_driver(2)
H = DenseHamiltonian([(s)->1-s, (s)->s], [-Hd, Hp], unit=:ħ)


plot(H, range(0,1,length=100), 4, linewidth=2)
xlabel!("s")
ylabel!(L"P(s)")


tf = 20
ρ0 = (σi+σx)⊗(σi+σx)/4
coupling = ConstantCouplings([σz⊗σi, σi⊗σz], unit=:ħ)
annealing = Annealing(H, ρ0, bath=bath, coupling=coupling)
close_sol = solve_von_neumann(annealing, tf, alg = Tsit5(), abstol=1e-6, reltol=1e-6);


plot(close_sol, H, 1, range(0,tf,length=100), linewidth=2)
xlabel!("t")
ylabel!(L"P_G(s)")


t_axis = range(0,tf,length=100)
p_computational_basis = [real(diag(close_sol(s))) for s in t_axis]
p_computational_basis = hcat(p_computational_basis...)
plot(t_axis, p_computational_basis', linewidth=2, label=[L"\rho_{00}" L"\rho_{11}" L"\rho_{22}" L"\rho_{33}"])
xlabel!("t")
ylabel!(L"\rho")


tf = 20
U = solve_unitary(annealing, tf, alg = Tsit5(), abstol=1e-7, reltol=1e-7);
redfield_sol = solve_redfield(annealing, tf, U, alg = Tsit5(), abstol=1e-7, reltol=1e-7);


t_axis = range(0,tf,length=100)
p_computational_basis = [real(diag(redfield_sol(s))) for s in t_axis]
p_computational_basis = hcat(p_computational_basis...)
plot(t_axis, p_computational_basis', linewidth=2, label=[L"\rho_{00}" L"\rho_{11}" L"\rho_{22}" L"\rho_{33}"])
xlabel!("t")
ylabel!(L"\rho")


redfield_sol = solve_redfield(annealing, tf, U, alg = Tsit5(), abstol=1e-7, reltol=1e-7, callback=PositivityCheckCallback())

