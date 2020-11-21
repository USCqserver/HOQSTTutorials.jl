function build_adiabatic_frame_evolution(hp, stage; polaron=false)
    ilvl = init_lvl(hp)
    gap_lvl = ilvl == 1 ? 1 : 2
    s_axis = range(0, 1, length=1001)
    H = build_H(hp, stage, polaron=polaron)

    # calculate the position of crossing
    sm, = min_gap(H, gap_lvl, gap_lvl + 1)

    dH = build_dH(hp, stage, polaron=polaron)
    coupling = build_couplings(polaron=polaron)
    if stage == 1
        inds = sm < 0.5 ? findlast((x) -> x < sm, s_axis) : nothing
        proj = project_to_lowlevel(H, s_axis, coupling, dH, lvl=5, direction=:backward)
        adiabatic_H = build_adiabatic_H(proj)
        couplings = interp_couplings(proj, inds, sm, polaron, stage)
    elseif stage == 2
        w, v = eigen_decomp(H, 0.0, lvl=5)
        adiabatic_H = AdiabaticFrameHamiltonian((s) -> w[1:4], nothing)
        couplings = rotate_couplings(coupling, v[:, 1:4], polaron)
    elseif stage == 3
        inds = sm > 0.5 ? findlast((x) -> x < sm, s_axis) : nothing
        proj = project_to_lowlevel(H, s_axis, coupling, dH, lvl=5)
        adiabatic_H = build_adiabatic_H(proj)
        couplings = interp_couplings(proj, inds, sm, polaron, stage)
    else
        throw(ArgumentError("Stage $stage not defined."))
    end
    cb = build_diabatic_callback(sm, gap_lvl, stage)
    adiabatic_H, couplings, cb
end

function build_diabatic_callback(sm, swap_lvl, stage)
    if stage == 1
        if sm < 0.5
            op = swap_operator(swap_lvl)
            cb = InstPulseCallback([2 * τ1 * sm], (c, i) -> c .= op * c * op')
        else
            cb = nothing
        end
    elseif stage == 2
        cb = nothing
    elseif stage == 3
        if sm > 0.5
            op = swap_operator(swap_lvl)
            cb = InstPulseCallback([2 * τ1 * sm], (c, i) -> c .= op * c * op')
        else
            cb = nothing
        end
    else
        throw(ArgumentError("Stage $stage not defined."))
    end
    cb
end

function build_adiabatic_H(proj)
    D = construct_interpolations(proj.s, hcat(proj.ev...)[1:4, :], order=2)
    Dfun(s) = [D(i, s) for i in 1:4]
    AdiabaticFrameHamiltonian(Dfun, nothing)
end

function min_gap(H, i=1, j=2)
    gap = function (s)
        w, = eigen_decomp(H, s, lvl=4)
        w[j] - w[i]
    end
    res1 = optimize(gap, 0, 0.5)
    res2 = optimize(gap, 0.5, 1)
    v, ind = findmin([Optim.minimum(res1), Optim.minimum(res2)])
    Optim.minimizer([res1, res2][ind]), v
end

function swap_operator(swap_lvl)
    a1 = zeros(4)
    a2 = zeros(4)
    a1[swap_lvl] = 1
    a2[swap_lvl + 1] = 1
    op = a1 * a2' + a2 * a1'
    resi = [i for i = 1:4 if !(i in [swap_lvl, swap_lvl + 1])]
    for i in resi
        op[i, i] = 1
    end
    op
end