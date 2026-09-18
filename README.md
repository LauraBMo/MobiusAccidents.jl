# MobiusAccidents

Enumerate, classify, and render the **accidental** Möbius transformations of the
d-th roots of unity — the *sporadic semi-dihedral* maps of Siliciano/Rank2Forms.

Formally, an **accident** is a Möbius map `ψ ∈ PGL₂(ℂ)` that sends **≥ 4 but not
all** of the d-th roots of unity `Θ_d = {e^{2πik/d}}` to roots of unity, while not
lying in the dihedral group `D_d`. Every accident preserves the unit circle
setwise, so the invariant circle satisfies **C = ψ(C)** — that circle, and the
permutation ψ induces on the roots sitting on it, is what this package is about.

It is the accident-specific layer on top of two generic siblings:

| Package | Role |
|---|---|
| [`MobiusSphere`](https://github.com/LauraBMo/MobiusSphere) | maths: PGL₂, stereographic projection, `Mobius_to_rigid_sitting` |
| [`MobiusSpherePlots`](https://github.com/LauraBMo/MobiusSpherePlots) | generic POV-Ray render of any sphere motion `(v, θ, t)` |
| **`MobiusAccidents`** | **enumerate accidents · bridge each to a rigid motion · overlay roots + images** |

It changes nothing about *how* a sphere is rendered — it feeds the generic renderer
the right motion and the root-of-unity floor overlays. Core fixes to the render
(coordinate conventions, chirality) are inherited from the sibling, not forked.

## Install

Developed alongside the siblings in the shared `@MobiusSuite` environment:

```julia
using Pkg
Pkg.develop(path="~/.julia/dev/MobiusAccidents")
```

Needs `povray` (≥ 3.7) and `ffmpeg` on `PATH` to render (inherited from
`MobiusSpherePlots`); classification and the bridge need neither.

## Usage

```julia
using MobiusAccidents

accident_table()                 # census vs the published counts (all ✓)

accs = classify(5)               # Vector{QuasiDihedral}, max-overlap first

# Render the d=5 maximum-overlap accident (root dots + ψ-image rings) to /tmp:
res = render_accident(5; quality=:medium)
res.path                         # the output GIF/MP4
res.maps, res.motions            # what was rendered (vectors)

# Render several at once, each a clip, concatenated into one video:
render_accident(5; which=:all)                    # every accident (+ the dihedral rep)
render_dihedral(5)                                # all 10 symmetries of D_5, in sequence
render_dihedral(5; which=[2, 3])                  # just a couple

# The raw rigid motion of any map (needs `using MobiusSphere`):
using MobiusSphere
Mobius_to_rot_angle_sitting(mobius_map(accs[1]))  # (; v, θ, t, imag_error)
```

`render_accident(d; which, overlays, output, …)` and its sibling `render_dihedral(d; …)`
select one or more maps — `which` is `:max`/`:all`/`:dihedral`, an index, a vector of
indices, or explicit `QuasiDihedral`s — render each with the root overlays, and
concatenate multiple clips into one video. The rest forwards to `render_mobius_animation`.

### Overlays

The overlays trace the **action of ψ on the overlapping roots** — the `k` roots the
accident carries to roots (only those, not all `d`). Each gets a distinct colour.
The **source** is a large flat disc; its **image** is a smaller disc of the same
colour that is **clock-animated** — it starts on the source, rides the caustic (the
moving rainbow/grid) through the whole motion, and lands at `ψ(ω^k)`. So you watch
each root travel to its image along the deforming pattern. Discs are flat, so a
smaller one stays visible on top of a larger one. Restyle via `marker_kwargs`
(`source_r`, `image_r`, `colors`, …), forwarded to `accident_overlay_sdl`.
