# MobiusSphereAccidentals

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
| [`MobiusSphereVisual`](https://github.com/LauraBMo/MobiusSphereVisual) | generic POV-Ray render of any sphere motion `(v, θ, t)` |
| **`MobiusSphereAccidentals`** | **enumerate accidents · bridge each to a rigid motion · overlay roots + images** |

It changes nothing about *how* a sphere is rendered — it feeds the generic renderer
the right motion and the root-of-unity floor overlays. Core fixes to the render
(coordinate conventions, chirality) are inherited from the sibling, not forked.

## Install

Developed alongside the siblings in the shared `@MobiusSuite` environment:

```julia
using Pkg
Pkg.develop(path="~/.julia/dev/MobiusSphereAccidentals")
```

Needs `povray` (≥ 3.7) and `ffmpeg` on `PATH` to render (inherited from
`MobiusSphereVisual`); classification and the bridge need neither.

## Usage

```julia
using MobiusSphereAccidentals

accident_table()                 # census vs the published counts (all ✓)

accs = classify(5)               # Vector{Accident}, max-overlap first
m = accident_to_rigid(accs[1])   # (; v, θ, t, imag_error) sphere rigid motion

# Render the d=5 maximum-overlap accident (root dots + ψ-image rings) to /tmp:
res = render_accident(5; quality=:medium)
res.path                         # the output GIF/MP4
res.accident, res.motion         # what was rendered
```

`render_accident(d; which, overlays, output, fps, nframes, resolution, quality, …)`
selects the representative (`:maxk`, an index, or a canonical exponent set), toggles
the overlays, and forwards the rest to `render_mobius_animation`.

### Overlays

The overlays trace the **action of ψ on the overlapping roots** — the `k` roots the
accident carries to roots (only those, not all `d`). Each gets a distinct colour;
its **source** is a large dot and its **image** `ψ(ω^k)` a smaller dot of the same
colour just above. A fixed root reads as one dot; a moved root shows a large dot
with a same-coloured small dot at its destination — the permutation made visible.
All sit on the invariant circle `|z|=1`. Restyle via `marker_kwargs`
(`source_size`, `image_size`, `colors`, …), forwarded to `accident_markers`.
