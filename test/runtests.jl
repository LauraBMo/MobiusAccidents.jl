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

    @testset "classify is sorted with :maxk first" begin
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
            @test ψ(Θ[s+1]) ≈ Θ[t+1] atol=1e-8
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

    @testset "accident_markers: colored source→image pairs on overlapping roots" begin
        a = classify(5)[1]                             # k = 4
        ms = accident_markers(a)
        @test length(ms) == 2 * a.k                    # a source + an image per overlap root
        @test all(m -> m.kind === :dot, ms)
        # every marker sits on the invariant circle of radius 1 (x² + z² = 1)
        for m in ms
            @test hypot(m.pos[1], m.pos[3]) ≈ 1.0 atol=1e-8
            @test m.pos[2] > 0                          # lifted above the floor
        end
        sizes = unique(round.([m.size for m in ms]; digits=6))
        @test length(sizes) == 2                        # one large (source), one small (image)
        # each colour appears exactly twice: its large source dot and small image dot
        cols = [m.color for m in ms]
        @test all(c -> count(==(c), cols) == 2, unique(cols))
        @test length(unique(cols)) == a.k               # k distinct colours
    end

    @testset "_select forms agree" begin
        accs = classify(8)
        a_max = MobiusSphereAccidentals._select(accs, :maxk)
        @test a_max.k == maximum(a.k for a in accs)
        @test MobiusSphereAccidentals._select(accs, 1) === accs[1]
        # canonical-set selector round-trips
        a1 = accs[1]
        @test MobiusSphereAccidentals._select(accs, a1.srcexp).srcexp == a1.srcexp
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
