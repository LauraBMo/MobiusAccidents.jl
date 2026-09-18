using MobiusAccidents
using MobiusSphere            # Mobius_to_rot_angle_sitting
using LinearAlgebra
using Test

@testset "MobiusAccidents.jl" begin

    @testset "classification census matches the published table" begin
        # (#accidents mod D_d, max overlap k) — the paper's counts.
        expected = Dict(5=>(1,4), 6=>(2,4), 7=>(4,4), 8=>(7,6), 9=>(10,4),
                        10=>(15,6), 11=>(20,4), 12=>(28,8), 13=>(35,4))
        for (d, (cnt, mk)) in sort(collect(expected))
            accs = classify(d)
            quasi = filter(!isdihedral, accs)             # the accidents (k < d)
            @test length(quasi) == cnt
            @test maximum(a.k for a in quasi) == mk
            @test all(a -> 4 <= a.k < d, quasi)           # accidents: ≥4 but not all
            @test count(isdihedral, accs) == 1            # plus one genuine D_d representative
        end
    end

    @testset "dihedral case admitted (k = d, ψ ∈ D_d)" begin
        for d in (5, 7, 12)
            accs = classify(d)
            dih = filter(isdihedral, accs)
            @test length(dih) == 1                        # orbit-dedup ⇒ a single rep
            @test dih[1].k == d == dih[1].d               # carries all d roots
            @test dih[1].srcexp == collect(0:d-1)         # overlap is the whole root set
            @test isdihedral(accs[end]) && !isdihedral(accs[1])  # sorted last, accidents first
        end
    end

    @testset "classify: accidents first (:max), dihedral last" begin
        accs = classify(12)
        @test !isdihedral(accs[1])                                   # an accident leads
        @test accs[1].k == maximum(a.k for a in accs if !isdihedral(a))
        @test isdihedral(accs[end])                                  # dihedral sorts last
        @test issorted(accs; by = a -> (isdihedral(a), -a.k, a.srcexp))
    end

    @testset "mobius_map realizes the source→target triple" begin
        a = classify(5)[1]
        ψ = mobius_map(a)
        Θ = roots_of_unity(5)
        # ψ sends ω^src → ω^tgt for the defining triple
        for (s, t) in zip(a.src, a.tgt)
            @test ψ(Θ[s]) ≈ Θ[t] atol=1e-8       # a.src/a.tgt are 1-based indices into Θ
        end
        # and carries exactly `k` roots to roots
        hits = count(e -> any(r -> abs(ψ(Θ[e+1]) - r) < 1e-8, Θ), 0:4)
        @test hits == a.k
    end

    @testset "Mobius_to_rot_angle_sitting: real motion + d=5 regression" begin
        a = classify(5)[1]                 # THE d=5 accident (unique mod D_5)
        m = Mobius_to_rot_angle_sitting(mobius_map(a))
        @test norm(m.v) ≈ 1.0 atol=1e-10   # unit axis
        # Pinned values (verified 2026-09-01 against the standalone bridge):
        @test m.v ≈ [0.7499818036699579, 0.5448936755963917, 0.3749909018349791] atol=1e-8
        @test rad2deg(m.θ) ≈ 81.8161275081839 atol=1e-6
        @test m.t ≈ [0.9756836610416143, -0.7088756736266997, -0.9213106741667365] atol=1e-8
    end

    @testset "dihedral_maps enumerates D_d as renderable motions" begin
        for d in (5, 6, 7)
            dm = dihedral_maps(d)
            @test length(dm) == 2d                          # |D_d| = 2d
            @test all(isdihedral, dm)                       # every one has k = d
            @test all(m -> m.srcexp == collect(0:d-1), dm)  # carries all roots
            @test dm[1].tgtexp == collect(0:d-1)            # identity r⁰ is first
            # a non-identity element realizes a genuine real sphere motion
            @test norm(Mobius_to_rot_angle_sitting(mobius_map(dm[2])).v) ≈ 1.0 atol=1e-10
        end
    end

    @testset "_cap_point: inverse stereo lands the rest shadow on z" begin
        # cap point rides on the unit sphere centred at <0,1,0>, on the lower cap,
        # and its shadow from the rest pole <0,2,0> projects back to z.
        for z in (cis(0.0), cis(1.3), cis(2.7), roots_of_unity(5)...)
            P = MobiusAccidents._cap_point(z)
            @test hypot(P[1], P[2] - 1, P[3]) ≈ 1.0 atol=1e-10   # on the sphere
            @test P[2] < 1                                        # lower cap
            λ = -2 / (P[2] - 2)                                   # project from <0,2,0> to y=0
            @test λ * P[1] ≈ real(z) atol=1e-10
            @test λ * P[3] ≈ imag(z) atol=1e-10
        end
    end

    @testset "accident_overlay_sdl: all d roots, source + animated image discs" begin
        a = classify(5)[1]                              # k = 4, d = 5
        sdl = accident_overlay_sdl(a)
        @test sdl isa AbstractString
        # every root gets a source + an image disc — all d of them (the k overlapping
        # roots in distinct hues, the rest as black source / white image discs)
        @test count("cylinder {", sdl) == 2 * a.d        # a source + an image disc per root
        @test count("vaxis_rotate", sdl) == a.d          # each image disc is clock-animated
        @test occursin("color rgb <0.7, 0.7, 0.7>", sdl) # non-overlapping source disc: grey
        @test occursin("color rgb <1.0, 1.0, 1.0>", sdl) # non-overlapping image disc: white
        @test occursin("PoleNow", sdl)                   # projected from the moving pole
        @test occursin("select(clock", sdl)              # phase-aware clock motion
        # wrong-length colour override (colours apply to the k overlapping roots) is rejected
        @test_throws ArgumentError accident_overlay_sdl(a; colors=[(1.0,0.0,0.0)])
    end

    @testset "_select_maps forms" begin
        accs = classify(8)
        S = MobiusAccidents._select_maps
        a_max = only(S(accs, :max))
        @test !isdihedral(a_max)                                     # :max is an accident
        @test a_max.k == maximum(a.k for a in accs if !isdihedral(a))
        @test S(accs, 1) == [accs[1]]                                # Integer index
        @test S(accs, [1, 2]) == accs[1:2]                           # Vector of indices
        @test S(accs, :all) === accs                                 # :all is the whole pool
        @test all(isdihedral, S(accs, :dihedral))                    # :dihedral filters
        @test S(accs, accs[1]) == [accs[1]]                          # explicit map
    end
end

# ── End-to-end render (requires povray + ffmpeg) ──────────────────────────────
if success(`which povray`) && success(`which ffmpeg`)
    @testset "render_accident end-to-end (:draft, tiny)" begin
        out = tempname() * ".gif"
        res = render_accident(5; quality=:draft, nframes=3, fps=3,
                              resolution=(320,180), output=out)
        @test isfile(res.path)
        @test filesize(res.path) > 1000
        @test length(res.maps) == 1
        @test res.maps[1].k == 4
        @test norm(res.motions[1].v) ≈ 1.0 atol=1e-10
        rm(res.path; force=true)
        println("  ✓  render_accident(5) → $(res.path)")
    end

    @testset "render_dihedral multi-render (two D_5 symmetries, concat)" begin
        out = tempname() * ".gif"
        res = render_dihedral(5; which=[2, 3], overlays=false, quality=:draft,
                              nframes=2, fps=3, resolution=(160,90), output=out)
        @test isfile(res.path)
        @test filesize(res.path) > 1000
        @test length(res.maps) == 2 && all(isdihedral, res.maps)   # two dihedral clips joined
        rm(res.path; force=true)
        println("  ✓  render_dihedral(5; which=[2,3]) → $(res.path)")
    end
else
    @warn "Skipping render end-to-end tests — povray/ffmpeg not both in PATH"
end

println("\nAll done.")
