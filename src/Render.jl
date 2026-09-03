# Render.jl — the bridge from an accident's Möbius map ψ to a sphere rigid motion,
# the root-of-unity floor overlays, and the high-level `render_accident` entry point.

"""
    accident_to_rigid(a::QuasiDihedral) -> (; v, θ, t, imag_error)

Convert the accident's Möbius map `ψ` to the rigid motion of the **sitting** sphere
(rests on the floor, centre one radius up) that the render realizes: rotation axis
`v` (unit, z-up MobiusSphere convention), angle `θ` (radians), translation `t`.

Uses `MobiusSphere.Mobius_to_rigid_sitting` — NOT the origin-centred `Mobius_to_rigid`
— because the scene's floor draws the unit circle at radius 1 (the centred variant's
invariant circle sits at radius 2 and would deform the drawn circle). `imag_error`
reports how far the recovered motion strayed from real (should be ~0 for an accident).
"""
function accident_to_rigid(a::QuasiDihedral)
    ψ = mobius_map(a)
    Q, T = Mobius_to_rigid_sitting(ψ)
    imerr = max(maximum(abs, imag.(Q)), maximum(abs, imag.(collect(T))))
    v, θ = rotation_axis_angle(real.(Q))
    return (v = Float64.(real.(collect(v))), θ = Float64(real(θ)),
            t = Float64.(real.(collect(T))), imag_error = imerr)
end

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

_select(accs::Vector{QuasiDihedral}, which::Symbol) =
    which === :max ? argmax(a -> a.k, Iterators.filter(!isdihedral, accs)) :
    error("unknown selector :$which (use :max or an Integer index)")

_select(accs::Vector{QuasiDihedral}, i) = accs[i]

# function _select(accs::Vector{QuasiDihedral}, I::AbstractVector{<:Integer})
#     d = accs[1].d
#     target = _canonical(collect(Int, set), d)
#     i = findfirst(a -> _canonical(a.srcexp, a.d) == target, accs)
#     i === nothing && error("no accident with canonical source set $(collect(set)) for d=$d")
#     return accs[i]
# end

"""
    render_accident(d; which=:max, overlays=true, output="/tmp/accident_d\$d.gif",
                    fps=20, nframes=60, resolution=(640,360), quality=:medium,
                    keep_temp=false, marker_kwargs=(;), kwargs...)
        -> (; path, accident, motion)

Enumerate the accidents of the `d`-th roots of unity, pick one, convert it to a
sphere rigid motion, and render the animation with `MobiusSphereVisual`.

`which` selects the representative: `:max` (default — a maximum-overlap accident) or
an `Integer` index into `classify(d)`. `overlays=true` draws the coloured source discs and their
clock-animated ψ-image discs (`marker_kwargs` forwards styling to
[`accident_overlay_sdl`](@ref)). Extra `kwargs...` pass straight through to
`render_mobius_animation` (e.g. `sampling`).

Returns the output `path` together with the chosen `accident` and its `motion`
`(; v, θ, t, imag_error)`, so the numbers are available for inspection.
"""
function render_accident(d::Integer;
        which = :max,
        overlays::Bool = true,
        output::AbstractString = "/tmp/accident_d$(d).gif",
        fps::Int = 20, nframes::Int = 60,
        resolution::Tuple{Int,Int} = (640, 360),
        quality::Symbol = :medium, keep_temp::Bool = false,
        marker_kwargs = (;), kwargs...)
    accs = classify(d)
    isempty(accs) && error("no accidents for d=$d")
    a = _select(accs, which)
    # if 'a' is a vector accs[[1,2,3]] : create files for 1,2 and 3 and concat them into output.
    m = accident_to_rigid(a)
    markers = overlays ? [(; kind = :raw, sdl = accident_overlay_sdl(a; marker_kwargs...))] : nothing
    path = render_mobius_animation(m.v, m.θ, m.t;
        output = String(output), fps = fps, nframes = nframes,
        resolution = resolution, quality = quality, keep_temp = keep_temp,
        markers = markers, kwargs...)
    return (; path, accident = a, motion = m)
end
