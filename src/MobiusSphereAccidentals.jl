"""
    MobiusSphereAccidentals

Enumerate, classify, and render the **accidental** Möbius transformations of the
d-th roots of unity — the paper's *sporadic semi-dihedral* maps: `ψ ∈ PGL₂(ℂ)`
with `|ψ(Θ_d) ∩ Θ_d| ≥ 4` but `ψ ∉ D_d` (they send ≥4 but not all of the d-th
roots of unity to roots of unity, and preserve the unit circle `|z|=1` setwise, so
the invariant circle satisfies `C = ψ(C)`).

This is the accident-specific layer on top of the generic siblings:

- [`MobiusSphere`](https://github.com/LauraBMo/MobiusSphere) — the maths (PGL₂,
  stereographic projection, `Mobius_to_rigid_sitting`);
- [`MobiusSphereVisual`](https://github.com/LauraBMo/MobiusSphereVisual) — the
  generic POV-Ray render of any `(v, θ, t)` sphere motion.

It adds nothing to how a sphere is rendered; it enumerates accidents, bridges each
to a sphere rigid motion, and feeds the generic renderer floor overlays marking the
roots of unity and their images.

# Quick start
```julia
using MobiusSphereAccidentals
accident_table()                 # census vs the published counts
render_accident(5)               # render the d=5 max-overlap accident to /tmp
```
"""
module MobiusSphereAccidentals

using LinearAlgebra
using MobiusSphere            # Möbius, Mobius_to_rigid_sitting, rotation_axis_angle
using MobiusSphereVisual      # render_mobius_animation

export Accident,
       roots_of_unity, classify, accident_table, mobius_map,
       accident_to_rigid, accident_markers, render_accident

include("Classify.jl")
include("Render.jl")

end # module
