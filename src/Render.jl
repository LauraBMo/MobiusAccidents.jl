# Render.jl — the bridge from a QuasiDihedral's Möbius map ψ to a sphere rigid motion
# (via MobiusSphere.Mobius_to_rot_angle_sitting), the root-of-unity floor overlays, and
# the high-level `render_accident` / `render_dihedral` entry points.
#
# We use the **sitting** decomposition (sphere resting on the floor, centre one radius
# up) — NOT the origin-centred one — because the scene's floor draws the unit circle at
# radius 1; the centred variant's invariant circle sits at radius 2 and would deform it.
# `Mobius_to_rot_angle_sitting` returns `(; v, θ, t)`: axis (unit, z-up), angle (radians),
# translation. There is no `imag_error` term — every Möbius map is realized by a genuine
# real rigid motion, so the deviation from real is identically zero and carries no signal.
_motion(a::QuasiDihedral) = Mobius_to_rot_angle_sitting(mobius_map(a))

# Evenly-spaced fully-saturated hue → RGB (the standard sat=1, val=1 hue ramp).
function _hue_rgb(h::Real)
    r = clamp(abs(6h - 3) - 1, 0, 1)
    g = clamp(2 - abs(6h - 2), 0, 1)
    b = clamp(2 - abs(6h - 4), 0, 1)
    return Float64[r, g, b]
end

# Inverse stereographic (sitting sphere, centre <0,1,0>, radius 1): the world cap
# point whose shadow — cast from the rest pole <0,2,0> — lands on the floor at the
# complex point z. This is the point on the shell that "carries" z's colour; letting
# it ride the render's Motion and re-projecting from the moving pole tracks z's image
# through the caustic. For |z|=1 this is <0.8·Re z, 0.4, 0.8·Im z>.
function _cap_point(z::Number)
    u, v = real(z), imag(z)
    ρ2 = u*u + v*v
    return (4u / (ρ2 + 4), 1 + (ρ2 - 4) / (ρ2 + 4), 4v / (ρ2 + 4))
end
    
_disc_sdl(x, ylo, yhi, r, c) = string(
    "cylinder { <", x[1], ", ", ylo, ", ", x[2], ">, <", x[1], ", ", yhi, ", ", x[2], ">, ", r,
    " pigment { color rgb <", c[1], ", ", c[2], ", ", c[3],
    "> } finish { ambient 1 diffuse 0 } no_shadow }")

"""
    accident_overlay_sdl(a::QuasiDihedral; source_r, image_r, source_h, image_h, colors)
        -> String

POV-Ray SDL (for a `:raw` marker) drawing the **action of ψ on the overlapping
roots** — the `k` roots the accident carries to roots (only those, not all `d`).
Each gets a distinct colour. For root `i`:

- a **large flat source disc** at `ω^srcexp[i]` — static;
- a **smaller flat image disc**, same colour, that is **clock-animated**: it starts
  on the source disc, rides the caustic (the moving rainbow/grid) through the whole
  motion, and lands at the image `ψ(ω^srcexp[i]) = ω^tgtexp[i]`.

The image disc's position is the shadow of root `i`'s cap point ([`_cap_point`](@ref))
after the scene's `Motion`, projected from the moving projector light `PoleNow` — so
it coincides with the caustic point of that root at every frame. Discs are flat (thin
cylinders) so a smaller disc stays visible on top of a larger one. `source_h`/`image_h`
are `(low, high)` disc heights (image sits above source); default colours are
evenly-spaced hues (override with `colors`, a length-`k` vector of RGB triples).
"""
function accident_overlay_sdl(a::QuasiDihedral;
        source_r::Real = 0.11, image_r::Real = 0.06,
        source_h = (0.008, 0.028), image_h = (0.032, 0.055),
        colors = nothing)
    d = a.d
    k = a.k
    Θ = roots_of_unity(d)
    cols = colors === nothing ? [_hue_rgb((i - 1) / k) for i in 1:k] :
                                [collect(Float64.(c)) for c in colors]
    length(cols) == k || throw(ArgumentError("colors must have length k=$k, got $(length(cols))"))
    print_src(IO, r, c, rad = source_r) =
        println(IO,
                _disc_sdl((real(r), imag(r)), source_h[1], source_h[2], rad, c))
    # Replicate the template's clock motion on the cap point, then project it from
    # the moving pole onto the floor (y=0): Fp is exactly this root's caustic point.
    print_image(IO, P, c) = 
        println(IO, """
        #declare MkAng = select(clock - 0.5, clock * 2 * Th, Th);
        #declare Pm = vaxis_rotate(<$(P[1]), $(P[2]), $(P[3])> - SphC0, Vax, MkAng) + SphC0;
        #declare Pm = Pm + Tv * select(clock - 0.5, 0, (clock - 0.5) * 2);
        #declare Lm = -PoleNow.y / (Pm.y - PoleNow.y);
        #declare Fp = PoleNow + Lm * (Pm - PoleNow);
        cylinder { <Fp.x, $(image_h[1]), Fp.z>, <Fp.x, $(image_h[2]), Fp.z>, $(image_r)
          pigment { color rgb <$(c[1]), $(c[2]), $(c[3])> } finish { ambient 1 diffuse 0 } no_shadow }""")
    io = IOBuffer()
    println(io, "// --- accident overlay: source discs (static) + image discs (clock-animated) ---")
    for i in 1:k                     
        r = Θ[a.srcexp[i] + 1]
        P = _cap_point(r)
        c = cols[i]
        # static large source disc
        print_src(io, r, c)
        # moving small image disc
        print_image(io, P, c)
    end
    for i in 0:(d-1)
        if i ∉ a.srcexp
            r = Θ[i + 1]
            P = _cap_point(r)
            # static source disc: light grey (not black, so it stays distinct from the
            # black grid it may cross) and enlarged for visibility
            print_src(io, r, [0.7,0.7,0.7], source_r * 1.6)
            # moving small white image disc
            print_image(io, P, [1.0,1.0,1.0])
        end
    end

    return String(take!(io))
end

# Select a sub-list of maps to render, from a pool (`classify(d)` or `dihedral_maps(d)`).
# Always returns a `Vector{QuasiDihedral}`; the render entry points render each and
# concatenate. `which` may be a Symbol (:max — the maximum-overlap accident; :all;
# :dihedral), an Integer index, a Vector of indices, a single QuasiDihedral, or a Vector.
_select_maps(pool::Vector{QuasiDihedral}, which::Symbol) =
    which === :all      ? pool :
    which === :max      ? [argmax(a -> a.k, Iterators.filter(!isdihedral, pool))] :
    which === :dihedral ? filter(isdihedral, pool) :
    error("unknown selector :$which (use :max, :all, :dihedral, an Integer, or a vector)")
_select_maps(pool::Vector{QuasiDihedral}, i::Integer) = [pool[i]]
_select_maps(pool::Vector{QuasiDihedral}, I::AbstractVector{<:Integer}) = pool[I]
_select_maps(::Vector{QuasiDihedral}, q::QuasiDihedral) = [q]
_select_maps(::Vector{QuasiDihedral}, Q::AbstractVector{QuasiDihedral}) = collect(Q)

# Render a list of maps: each to its own clip (with its own overlay), then concat in
# order (see `concat_clips`). A single map renders straight to `output`. Per-map clips go
# to `<stem>_parts/part_NN.<ext>` (kept only with `keep_temp`). Returns the output `path`,
# the chosen `maps`, and their `motions` `(; v, θ, t)`.
function _render_maps(maps::Vector{QuasiDihedral};
        overlays::Bool = true,
        output::AbstractString = "/tmp/mobius.gif",
        fps::Int = 20, nframes::Int = 60,
        resolution::Tuple{Int,Int} = (640, 360),
        quality::Symbol = :medium, keep_temp::Bool = false,
        marker_kwargs = (;), kwargs...)
    isempty(maps) && error("no maps to render")
    motions = [_motion(a) for a in maps]
    overlay(a) = overlays ? [(; kind = :raw, sdl = accident_overlay_sdl(a; marker_kwargs...))] : nothing
    render1(a, m, out) = render_mobius_animation(m.v, m.θ, m.t;
        output = out, fps = fps, nframes = nframes, resolution = resolution,
        quality = quality, keep_temp = keep_temp, markers = overlay(a), kwargs...)

    if length(maps) == 1
        path = render1(maps[1], motions[1], String(output))
        return (; path, maps, motions)
    end

    stem, ext = splitext(String(output))
    partdir = stem * "_parts"
    mkpath(partdir)
    parts = String[]
    for (i, a) in enumerate(maps)
        req = joinpath(partdir, "part_" * lpad(i, 2, '0') * ext)
        push!(parts, render1(a, motions[i], req))   # actual path (may differ on format fallback)
    end
    # parts may have fallen back to another format (e.g. mp4→gif); match output to them.
    path = concat_clips(parts, stem * splitext(parts[1])[2])
    keep_temp || rm(partdir; recursive = true, force = true)
    return (; path, maps, motions)
end

"""
    render_accident(d; which=:max, overlays=true, output="/tmp/accident_d\$d.gif",
                    fps=20, nframes=60, resolution=(640,360), quality=:medium,
                    keep_temp=false, marker_kwargs=(;), kwargs...)
        -> (; path, maps, motions)

Enumerate the accidents of the `d`-th roots of unity ([`classify`](@ref)), pick one or
several, convert each to a sitting-sphere rigid motion, and render with
`MobiusSpherePlots`. `which` selects from `classify(d)`: `:max` (default — the
maximum-overlap accident), `:all`, `:dihedral`, an `Integer` index, a `Vector{Int}` of
indices, or explicit `QuasiDihedral`(s). A multi-map selection renders each to its own
clip and concatenates them ([`concat_clips`](@ref)).

`overlays=true` draws the coloured root discs + their clock-animated ψ-image discs
(`marker_kwargs` forwards styling to [`accident_overlay_sdl`](@ref)). Extra `kwargs...`
pass through to `render_mobius_animation` (e.g. `sampling`, `hold`). Returns
`(; path, maps, motions)`.
"""
function render_accident(d::Integer; which = :max,
        output::AbstractString = "/tmp/accident_d$(d).gif", kwargs...)
    accs = classify(d)
    isempty(accs) && error("no accidents for d=$d")
    return _render_maps(_select_maps(accs, which); output = output, kwargs...)
end

"""
    render_dihedral(d; which=:all, overlays=true, output="/tmp/dihedral_d\$d.gif", …)
        -> (; path, maps, motions)

Like [`render_accident`](@ref) but over the genuine dihedral symmetries
[`dihedral_maps(d)`](@ref) — the `2d` elements of `D_d` (rotations `z ↦ ω^k z`,
reflections `z ↦ ω^k/z`). `which` defaults to `:all` (every symmetry, concatenated;
the identity `k=0` is first and renders static). Same overlays/kwargs as
`render_accident`.
"""
render_dihedral(d::Integer; which = :all,
        output::AbstractString = "/tmp/dihedral_d$(d).gif", kwargs...) =
    _render_maps(_select_maps(dihedral_maps(d), which); output = output, kwargs...)
