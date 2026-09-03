using MobiusSphereAccidentals
using LinearAlgebra
using Test

@testset "MobiusSphereAccidentals.jl" begin

    @testset "classification census matches the published table" begin
        # (#accidents mod D_d, max overlap k) — the paper's counts.
        expected = Dict(5=>(1,4), 6=>(2,4), 7=>(4,4), 8=>(7,6), 9=>(10,4),
                        10=>(15,6), 11=>(20,4), 12=>(28,8), 13=>(35,4))
        for (d, (cnt, mk)) in sort(collect(expected))
            accs = classify(d)
            @test length(accs) == cnt
            @test maximum(a.k for a in accs) == mk
            @test all(a -> 4 <= a.k < d, accs)      # accidents: ≥4 but not all
        end
    end

    @testset "classify is sorted with :max first" begin
        accs = classify(12)
        @test accs[1].k == maximum(a.k for a in accs)   # max-overlap leads
        @test issorted(accs; by = a -> (-a.k, a.srcexp))
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

    @testset "accident_to_rigid: real motion + d=5 regression" begin
        a = classify(5)[1]                 # THE d=5 accident (unique mod D_5)
        m = accident_to_rigid(a)
        @test m.imag_error < 1e-8          # ψ gives a genuine real rigid motion
        @test norm(m.v) ≈ 1.0 atol=1e-10   # unit axis
        # Pinned values (verified 2026-09-01 against the standalone bridge):
        @test m.v ≈ [0.7499818036699579, 0.5448936755963917, 0.3749909018349791] atol=1e-8
        @test rad2deg(m.θ) ≈ 81.8161275081839 atol=1e-6
        @test m.t ≈ [0.9756836610416143, -0.7088756736266997, -0.9213106741667365] atol=1e-8
    end

    @testset "_cap_point: inverse stereo lands the rest shadow on z" begin
        # cap point rides on the unit sphere centred at <0,1,0>, on the lower cap,
        # and its shadow from the rest pole <0,2,0> projects back to z.
        for z in (cis(0.0), cis(1.3), cis(2.7), roots_of_unity(5)...)
            P = MobiusSphereAccidentals._cap_point(z)
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

    @testset "_select forms agree" begin
        accs = classify(8)
        a_max = MobiusSphereAccidentals._select(accs, :max)
        @test a_max.k == maximum(a.k for a in accs)
        @test MobiusSphereAccidentals._select(accs, 1) === accs[1]
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
        @test res.accident.k == 4
        @test res.motion.imag_error < 1e-8
        rm(res.path; force=true)
        println("  ✓  render_accident(5) → $(res.path)")
    end
else
    @warn "Skipping render_accident end-to-end test — povray/ffmpeg not both in PATH"
end

println("\nAll done.")
