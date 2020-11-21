function build_couplings(; polaron=false)
    σ₊ = [0.0 1;0 0]
    σ₋ = [0.0 0;1 0]
    if polaron == false
        ConstantCouplings(["10ZII", "IZI", "IIZ"], unit=:ħ)
    else
        ConstantCouplings([σ₊ ⊗ σi ⊗ σi, σ₋ ⊗ σi ⊗ σi, σi ⊗ σz ⊗ σi, σi ⊗ σi ⊗ σz], unit=:ħ)
    end
end

function interp_couplings(proj, inds, sm, polaron, stage)
    if polaron == false
        ops = [[op[i][1:4, 1:4] for op in proj.op] for i in 1:3]
    else
        spfun, = build_s(stage)
        As = A.(spfun.(proj.s))
        ops1 = [As .* [op[i][1:4, 1:4] for op in proj.op] for i in 1:2]
        ops2 = [[op[i][1:4, 1:4] for op in proj.op] for i in 3:4]
        ops = vcat(ops1, ops2)
    end
    ops = [interp_single_coupling(proj.s, x, inds, sm) for x in ops]
    if polaron == false
        [CustomCouplings(ops, unit=:ħ), CustomCouplings([ops[1]], unit=:ħ)]
    else
        [CustomCouplings(ops[1:2], unit=:ħ), CustomCouplings(ops[3:4], unit=:ħ)]
    end
end

function interp_single_coupling(s, op, ::Nothing, ::Any)
    coupling = function (s, itp)
        res = Array{Float64,2}(undef, 4, 4)
        for i = 1:4
            for j = 1:4
                res[i, j] = itp(i, j, s)
            end
        end
        res
    end
    itp = construct_interpolations(s, op)
    (s) -> coupling(s, itp)
end

function interp_single_coupling(s_axis, op, inds, sm)
    coupling = function (s, itp)
        res = Array{Float64,2}(undef, 4, 4)
        for i = 1:4
            for j = 1:4
                res[i, j] = itp(i, j, s)
            end
        end
        res
    end
    itp1 = construct_interpolations(s_axis[1:inds], op[1:inds])
    itp2 = construct_interpolations(s_axis[inds + 1:end], op[inds + 1:end])
    function (s)
        if s <= sm
            coupling(s, itp1)
        else
            coupling(s, itp2)
        end
    end
end

function rotate_couplings(coupling, v, polaron)
    cmat = coupling(0)
    cmat = [real.(v' * m * v) for m in cmat]
    if polaron == false
        [ConstantCouplings(cmat, unit=:ħ), ConstantCouplings([cmat[1]], unit=:ħ)]
    else
        [ConstantCouplings(A(sp) * cmat[1:2], unit=:ħ), ConstantCouplings(cmat[3:4], unit=:ħ)]
    end
end