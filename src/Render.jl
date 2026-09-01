# Render.jl — the bridge from an accident's Möbius map ψ to a sphere rigid motion,
# the root-of-unity floor overlays, and the high-level `render_accident` entry point.

"""
    accident_to_rigid(a::Accident) -> (; v, θ, t, imag_error)

Convert the accident's Möbius map `ψ` to the rigid motion of the **sitting** sphere
(rests on the floor, centre one radius up) that the render realizes: rotation axis
`v` (unit, z-up MobiusSphere convention), angle `θ` (radians), translation `t`.

Uses `MobiusSphere.Mobius_to_rigid_sitting` — NOT the origin-centred `Mobius_to_rigid`
— because the scene's floor draws the unit circle at radius 1 (the centred variant's
invariant circle sits at radius 2 and would deform the drawn circle). `imag_error`
reports how far the recovered motion strayed from real (should be ~0 for an accident).
"""
function accident_to_rigid(a::Accident)
    ψ = mobius_map(a)
    Q, T = Mobius_to_rigid_sitting(ψ)
    imerr = max(maximum(abs, imag.(Q)), maximum(abs, imag.(collect(T))))
    v, θ = rotation_axis_angle(real.(Q))
    return (v = Float64.(real.(collect(v))), θ = Float64(real(θ)),
            t = Float64.(real.(collect(T))), imag_error = imerr)
end

# Floor ↔ complex plane: the caustic realizes ψ on the floor (y = 0 in POV world),
# with complex z = u + iv drawn at world <u, 0, v> and the unit circle at radius 1.
# `lift` nudges markers just above the floor so they don't z-fight the plane.
_floor_pos(z::Number; lift::Real = 0.02) = Float64[real(z), Float64(lift), imag(z)]

# Evenly-spaced fully-saturated hue → RGB (the standard sat=1, val=1 hue ramp).
function _hue_rgb(h::Real)
    r = clamp(abs(6h - 3) - 1, 0, 1)
    g = clamp(2 - abs(6h - 2), 0, 1)
    b = clamp(2 - abs(6h - 4), 0, 1)
    return Float64[r, g, b]
end

"""
    accident_markers(a::Accident; source_size, image_size, source_lift, image_lift,
                     colors) -> Vector{NamedTuple}

Floor markers for the `markers` keyword of `render_mobius_animation`, tracing the
**action of ψ on the overlapping roots** — the `k` roots the accident carries to
roots. Only those are marked (not all `d`). For each overlapping root:

- the **source** `ω^srcexp[i]` is a **large** dot;
- its **image** `ω^tgtexp[i] = ψ(ω^srcexp[i])` is a **smaller** dot in the **same
  colour**, lifted just above it.

So a fixed root reads as a single dot, and a moved root shows a large dot with a
same-coloured small dot at its destination — the permutation made visible. Each of
the `k` roots gets a distinct colour (evenly-spaced hues by default; override with
`colors`, a length-`k` vector of RGB triples). All markers sit on the invariant
circle `|z|=1` (radius 1 on the floor); they are static, and the caustic deforms
under them from identity to `ψ`.
"""
function accident_markers(a::Accident;
        source_size::Real = 0.09, image_size::Real = 0.05,
        source_lift::Real = 0.02, image_lift::Real = 0.06,
        colors = nothing)
    Θ = roots_of_unity(a.d)
    k = a.k
    cols = colors === nothing ? [_hue_rgb((i - 1) / k) for i in 1:k] :
                                [collect(Float64.(c)) for c in colors]
    length(cols) == k || throw(ArgumentError("colors must have length k=$k, got $(length(cols))"))
    ms = NamedTuple[]
    for i in 1:k       # large source dots first
        p = Θ[a.srcexp[i] + 1]
        push!(ms, (; kind = :dot, pos = _floor_pos(p; lift = source_lift),
                     color = cols[i], size = Float64(source_size)))
    end
    for i in 1:k       # smaller image dots on top (same colour as their source)
        p = Θ[a.tgtexp[i] + 1]
        push!(ms, (; kind = :dot, pos = _floor_pos(p; lift = image_lift),
                     color = cols[i], size = Float64(image_size)))
    end
    return ms
end

_select(accs::Vector{Accident}, which::Symbol) =
    which === :maxk ? argmax(a -> a.k, accs) :
    error("unknown selector :$which (use :maxk, an Integer index, or a Vector{Int} set)")

_select(accs::Vector{Accident}, i::Integer) = accs[i]

function _select(accs::Vector{Accident}, set::AbstractVector{<:Integer})
    d = accs[1].d
    target = _canonical(collect(Int, set), d)
    i = findfirst(a -> _canonical(a.srcexp, a.d) == target, accs)
    i === nothing && error("no accident with canonical source set $(collect(set)) for d=$d")
    return accs[i]
end

"""
    render_accident(d; which=:maxk, overlays=true, output="/tmp/accident_d\$d.gif",
                    fps=20, nframes=60, resolution=(640,360), quality=:medium,
                    keep_temp=false, marker_kwargs=(;), kwargs...)
        -> (; path, accident, motion)

Enumerate the accidents of the `d`-th roots of unity, pick one, convert it to a
sphere rigid motion, and render the animation with `MobiusSphereVisual`.

`which` selects the representative: `:maxk` (default — a maximum-overlap accident),
an `Integer` index into `classify(d)`, or a `Vector{Int}` giving a canonical
source-exponent set. `overlays=true` draws the root-of-unity dots and their
ψ-images (`marker_kwargs` forwards styling to [`accident_markers`](@ref)). Extra
`kwargs...` pass straight through to `render_mobius_animation` (e.g. `sampling`).

Returns the output `path` together with the chosen `accident` and its `motion`
`(; v, θ, t, imag_error)`, so the numbers are available for inspection.
"""
function render_accident(d::Integer;
        which = :maxk,
        overlays::Bool = true,
        output::AbstractString = "/tmp/accident_d$(d).gif",
        fps::Int = 20, nframes::Int = 60,
        resolution::Tuple{Int,Int} = (640, 360),
        quality::Symbol = :medium, keep_temp::Bool = false,
        marker_kwargs = (;), kwargs...)
    accs = classify(d)
    isempty(accs) && error("no accidents for d=$d")
    a = _select(accs, which)
    m = accident_to_rigid(a)
    markers = overlays ? accident_markers(a; marker_kwargs...) : nothing
    path = render_mobius_animation(m.v, m.θ, m.t;
        output = String(output), fps = fps, nframes = nframes,
        resolution = resolution, quality = quality, keep_temp = keep_temp,
        markers = markers, kwargs...)
    return (; path, accident = a, motion = m)
end
