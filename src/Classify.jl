# Classify.jl — Nemo-free enumeration and classification of the accidental Möbius
# transformations of the d-th roots of unity. Everything is plain ComplexF64;
# no cyclotomic field is needed. Ported from the (now retired) standalone script
# `accidental_mobius.jl`, which this package supersedes.

"""
    roots_of_unity(d) -> Vector{ComplexF64}

The d-th roots of unity `Θ_d = {exp(2πik/d) : k = 0,…,d-1}`.
"""
roots_of_unity(d::Integer) = [cispi(2 * k / d) for k in 0:d-1]

# Dihedral group D_d acting on exponent sets S ⊆ ℤ/dℤ, via DihedralGroups.jl. The
# group action `i^g` sends an exponent i by a rotation (i↦i+k) or a reflection
# (i↦k−i), both mod d; the orbit of S is its image under all 2d elements of D_d.
_Dd_orbit(S::Vector{Int}, d::Int) =
    [sort([i^g for i in S]) for g in dihedralgroup(d)]

_canonical(S, d) = minimum(_Dd_orbit(S, d))          # lexicographic orbit representative

## TODO TODO Use julia indeces instead of exponents,
## Θsrc = Θ[a.src] this calls work out of the box.
## TODO A QuasiDihedral may also represent a genuine Dihedral symmetry, when
## srcexp = Θ_d (k = d). Admit those too so we can plot them.
"""
    QuasiDihedral

One representative "accidental" (the paper's *sporadic semi-dihedral*) Möbius map
of the `d`-th roots of unity: a `ψ ∈ PGL₂(ℂ)` with `|ψ(Θ_d) ∩ Θ_d| ≥ 4` but
`ψ ∉ D_d`. Fields:

- `d`      — order of the root system
- `src`    — source **index** vector `[1, i+1, j+1]` into `Θ_d` defining `ψ` (1-based,
             so `Θ[a.src]` returns the defining source roots directly)
- `tgt`    — target index vector `[1, ip+1, jp+1]` into `Θ_d`, aligned with `src`
- `srcexp` — all source exponents `e` with `ψ(ω^e) ∈ Θ_d` (the overlap, sorted)
- `tgtexp` — their images' exponents, aligned with `srcexp`
- `k`      — overlap size `|ψ(Θ_d) ∩ Θ_d|` = `length(srcexp)`
"""
struct QuasiDihedral
    d::Int
    src::Vector{Int}
    tgt::Vector{Int}
    srcexp::Vector{Int}
    tgtexp::Vector{Int}
    k::Int
end

function Base.show(io::IO, a::QuasiDihedral)
    print(io, "QuasiDihedral(d=$(a.d), k=$(a.k), ",
          "src=$(a.srcexp) → tgt=$(sort(a.tgtexp)))")
end

"""
    mobius_map(a::QuasiDihedral)

The Möbius transformation `ψ` realizing the accident: it sends the source root
triple to the target root triple (and, being an accident, carries `k ≥ 4` of the
d-th roots of unity to roots of unity while preserving `|z|=1` setwise).
"""
function mobius_map(a::QuasiDihedral)
    Θ = roots_of_unity(a.d)
    Möbius(Θ[a.src], Θ[a.tgt])          # a.src / a.tgt are 1-based indices into Θ
end

"""
    classify(d; tol=1e-8) -> Vector{QuasiDihedral}

Enumerate every Möbius sending a canonical triplet `(1, ω^i, ω^j)` to
`(1, ω^ip, ω^jp)`, keep those whose overlap `|ψ(Θ_d) ∩ Θ_d|` is `≥ 4` and `< d`,
and dedup by the `D_d`-orbit of the overlap's source-exponent set. The result is
sorted by descending overlap `k`, then lexicographically by source set — so
`classify(d)[1]` (equivalently `which = :max`) is a maximum-overlap accident.
"""
function classify(d::Integer; tol = 1e-8)
    d = Int(d)
    Θ = roots_of_unity(d)
    onroot(w) = isfinite(w) ? findfirst(r -> abs(r - w) < tol, Θ) : nothing
    reps = Dict{Vector{Int},QuasiDihedral}()
    for i in 1:d-1, j in i+1:d-1, ip in 1:d-1, jp in ip+1:d-1
        ψ = Möbius([Θ[1], Θ[i+1], Θ[j+1]], [Θ[1], Θ[ip+1], Θ[jp+1]])
        srcexp, tgtexp = Int[], Int[]
        for e in 0:d-1
            k = onroot(ψ(Θ[e+1]))
            k === nothing && continue
            push!(srcexp, e); push!(tgtexp, k - 1)
        end
        K = length(srcexp)
        (4 <= K < d) || continue
        c = _canonical(srcexp, d)
        haskey(reps, c) ||
            (reps[c] = QuasiDihedral(d, [1, i+1, j+1], [1, ip+1, jp+1], sort(srcexp), tgtexp, K))
    end
    accs = collect(values(reps))
    sort!(accs; by = a -> (-a.k, a.srcexp))
    return accs
end

# ── Cross-check against the published table (count mod D_d, max overlap k) ──────────
const EXPECTED = Dict(5=>(1,4), 6=>(2,4), 7=>(4,4), 8=>(7,6), 9=>(10,4),
                      10=>(15,6), 11=>(20,4), 12=>(28,8), 13=>(35,4))

"""
    accident_table(ds = 5:13; io = stdout)

Print the accident census — number of accidents mod `D_d` and the maximum overlap
`k` — for each `d`, alongside the published expected values. A `✓` in the last
column means this enumeration agrees with the paper.
"""
function accident_table(ds = 5:13; io::IO = stdout)
    println(io, "  d | #accidents(modD) | max k | expected | ok")
    println(io, "----+------------------+-------+----------+----")
    for d in ds
        accs = classify(d)
        cnt = length(accs)
        mk = isempty(accs) ? 0 : maximum(a.k for a in accs)
        e = get(EXPECTED, d, nothing)
        ok = e === nothing ? "?" : ((cnt, mk) == e ? "✓" : "✗")
        estr = e === nothing ? "   --   " : lpad("($(e[1]),$(e[2]))", 8)
        println(io, lpad(d,3), " | ", lpad(cnt,16), " | ", lpad(mk,5), " | ", estr, " | ", ok)
    end
end
