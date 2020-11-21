const s_star = 0.339
const sp = 0.612
const τ1 = 10000
const vi = PauliVec[3][2] ⊗ PauliVec[3][2] ⊗ PauliVec[3][2]

"""
    build_s(stage)

Build the annealing parameter `s` based on `stage`. The stage can be 1, 2 or 3, as defined in HOQST paper.
"""
function build_s(stage)
    if stage == 1
        spfun = (s) -> s <= 0.5 ? 1.0 : 2 * (sp - 1) * s + 2 - sp
        sfun = (s) -> s <= 0.5 ? 2 * (s_star - 1) * s + 1 : s_star
    elseif stage == 2
        spfun = (s) -> sp
        sfun = (s) -> s_star
    elseif stage == 3
        spfun = (s) -> s <= 0.5 ? 2 * (1 - sp) * s + sp : 1.0
        sfun = (s) -> s <= 0.5 ? s_star : 2 * (1 - s_star) * s + 2 * s_star - 1
    else
        throw(ArgumentError("Stage $stage does not exist."))
    end
    spfun, sfun
end

"""
    build_ds(stage)

Build the derivative of the annealing parameter `s` based on `stage`.
"""
function build_ds(stage)
    if stage == 1
        dspfun = (s) -> s <= 0.5 ? 0.0 : 2 * (sp - 1)
        dsfun = (s) -> s <= 0.5 ? 2 * (s_star - 1) : 0.0
    elseif stage == 2
        dspfun = (s) -> 0.0
        dsfun = (s) -> 0.0
    elseif stage == 3
        dspfun = (s) -> s <= 0.5 ? 2 * (1 - sp) : 0.0
        dsfun = (s) -> s <= 0.5 ? 0.0 : 2 * (1 - s_star)
    else
        throw(ArgumentError("Stage $stage does not exists."))
    end
    dspfun, dsfun
end

function build_Hlist(hp; polaron=false)
    polaron == true ? 
    [-σi ⊗ σx ⊗ σi - σi ⊗ σi ⊗ σx, build_Hising(hp)] :
    [-σx ⊗ σi ⊗ σi, -σi ⊗ σx ⊗ σi - σi ⊗ σi ⊗ σx, build_Hising(hp)]
end

function build_Hising(hp)
    jp = -1.8
    js = -2.5
    hl = local_field_term([-hp, -jp, 0], [1, 2, 3], 3)
    h2 = two_local_term([jp, js], [[1, 2], [2, 3]], 3)
    hp = hl + h2
end

"""
    build_H(hp, stage; polaron=false)

Build the `stage` stage Hamiltonian given `hp`. `polaron` specifies whether the Hamiltonian is in polaron frame or not.
"""
function build_H(hp, stage; polaron=false)
    if stage == 1
        spfun, sfun = build_s(stage)
        Ap = (s) -> A(spfun(s))
        As = (s) -> A(sfun(s))
        Bs = (s) -> B(sfun(s))
        Hlist = build_Hlist(hp, polaron=polaron)
        polaron == false ? DenseHamiltonian([Ap, As, Bs], Hlist, unit=:ħ) : DenseHamiltonian([As, Bs], Hlist, unit=:ħ)
    elseif stage == 2
        hlist = build_Hlist(hp, polaron=polaron)
        hmat = polaron == false ? A(sp) * hlist[1] + A(s_star) * hlist[2] + B(s_star) * hlist[3] : A(s_star) * hlist[1] + B(s_star) * hlist[2]
        function build_user_eigen(u_cache)
            w, v = eigen(Hermitian(hmat))
            EIGS = function (H, t, lvl)
                w[1:lvl], v[:, 1:lvl]
            end
        end
        DenseHamiltonian([(s) -> 1.0], [hmat], unit=:ħ, EIGS=build_user_eigen)
    elseif stage == 3
        spfun, sfun = build_s(stage)
        Ap(s) = A(spfun(s))
        As(s) = A(sfun(s))
        Bs(s) = B(sfun(s))
        Hlist = build_Hlist(hp, polaron=polaron)
        polaron == false ? DenseHamiltonian([Ap, As, Bs], Hlist, unit=:ħ) : DenseHamiltonian([As, Bs], Hlist, unit=:ħ)
    else
        error("Stage $stage not defined.")
    end
end

"""
    build_dH(hp, stage; polaron=false)

Build the derivative of the `stage` stage Hamiltonian as a function.
"""
function build_dH(hp, stage; polaron=false)
    if stage in [1,2,3]
        spfun, sfun = build_s(stage)
        dspfun, dsfun = build_ds(stage)
        dAp(s) = dA(spfun(s)) * dspfun(s)
        dAs(s) = dA(sfun(s)) * dsfun(s)
        dBs(s) = dB(sfun(s)) * dsfun(s)
        Hlist = build_Hlist(hp, polaron=polaron)
        polaron == false ? (s) -> dAp(s) * Hlist[1] + dAs(s) * Hlist[2] + dBs(s) * Hlist[3] : (s) -> dAs(s) * Hlist[1] + dBs(s) * Hlist[2]
    else
        error("Stage $stage not defined.")
    end
end

"""
    init_lvl(hp)

Given `hp` value, check the initial energy level of the all one state.
"""
function init_lvl(hp)
    w, v = eigen(Hermitian(build_Hising(hp)))
    p = abs2.(v' * vi)
    ind = findfirst((x) -> x ≈ 1.0, p)
end

"""
    build_adiabatic_frame_u0(hp)

Build the adiabatic frame initial state given `hp` value.
"""
function build_adiabatic_frame_u0(hp)
    ilvl = init_lvl(hp)
    u0 = zeros(ComplexF64, 4)
    u0[ilvl] = 1.0
    u0
end

function hp_to_gap(hp)
    2 * B(s_star) * hp
end

function gap_to_hp(g)
    g / 2 / B(s_star)
end