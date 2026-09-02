# Joint NISAR + Sentinel-1 co-event decomposition (30-m reprocessed NISAR inputs)

- four geometries (NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23), decay-normalized to 18 July (factors 1.057, 1.032, 1.000, 1.105); NISAR from the 30-m reprocessed unwrapped fields

## Constrained decomposition (rift-perpendicular horizontal + vertical)

- two components per cell; the horizontal direction is fixed to the rift-perpendicular azimuth 63.0° taken from the fitted strike, so this variant is not independent of the source-model geometry
- cells solved: 3461 (two-look cells are exact solves; exactly three geometries: 1217; four geometries: 2020; total with at least three: 3237)
- flank horizontal probes (2-km disks): EAST +58.4 cm, WEST -55.8 cm ⇒ **opening divergence +114.2 cm** (uniform model: +90.8 cm)
- axial vertical probe: -13.0 cm (model +6.6 cm) — two-look cells only (both NISAR; no Sentinel-1 coverage survives at the axis)
- median formal sigmas: horizontal 1.0 cm, vertical 0.7 cm
- closure (weighted residual, σ units; redundant cells only): all median 0.58 (16% > 1σ, n=3237); within_25km median 0.82 (32% > 1σ, n=1109); within_15km median 1.17 (62% > 1σ, n=356); beyond_28km median 0.44 (7% > 1σ, n=1825)
- probe-disk closure (the disks that define the reported values; redundant cells only): EAST median 3.18 σ (n=6); WEST median 0.60 σ (n=4); AXIS no redundant cells

- zero-north variant (comparability with earlier products): divergence +102.9 cm, axis -7.3 cm
- no-decay-normalization variant: divergence +119.4 cm, axis -13.0 cm — the normalization is a model choice, and this is its size

## Free three-component decomposition (east, north, up; ≥3 looks)

- cells solved: 3237 (exactly three geometries: 1217; four geometries: 2020); no direction constraint
- median formal sigmas: east 0.9 cm, north 4.5 cm, up 0.8 cm — north is weakly observed by near-polar geometries
- rift-perpendicular flank probes projected from (E, N): EAST +59.2 cm, WEST -61.6 cm ⇒ divergence +120.8 cm
- axial vertical probe: decorrelated (three-look cells at best near the axis; 0 probe cells)
- probe-disk closure: EAST median 2.82 σ (n=5); WEST median 0.56 σ (n=3); AXIS no redundant cells

## Dike-frame components (rotated free three-component solve)

- frame from the fitted strike 153.0°: U∥ along strike (positive toward azimuth 333.0°, NNW — the near-north direction), U⊥ perpendicular-left (positive toward azimuth 63.0°, ENE), up; an exact rotation of the free solve, so coverage and closure are unchanged
- flank U⊥ probes: EAST +59.2 cm, WEST -61.6 cm ⇒ divergence +120.8 cm (identical to the projected values above)
- flank U∥ probes: EAST -5.3 cm, WEST +4.2 cm (uniform model: EAST -0.5 cm, WEST -0.2 cm)
- median formal sigmas: U⊥ 2.1 cm, U∥ 4.1 cm, up 0.8 cm — U∥ is nearly north–south, the direction the near-polar looks observe weakly, so it inherits most of the north uncertainty
- U∥ within 18 km: p01/p99 -39.6/+14.9 cm; 36% of cells exceed 2 σ∥
- strike sensitivity: with the southernmost 2005 segment strike 148.5° the flank probes become U⊥ EAST +58.6 cm, WEST -61.1 cm; U∥ EAST -9.9 cm, WEST +9.1 cm

## Scope

- The EAST-flank redundant cells close poorly (see the probe-disk
  lines): the looks disagree systematically on the east-lobe
  amplitude, consistent with the Sentinel-1-vs-NISAR amplitude
  ratios (0.88–0.92) in the S1 cross-checks. The joint horizontal
  value there is a cross-sensor compromise, not a closed
  measurement, and the divergence numbers inherit that systematic.
- The rift-perpendicular constraint is exact only if horizontal motion
  is opening-normal; along-strike motion maps into the closure, not
  into h. The free three-component solve has no such constraint but
  carries several-cm north sigmas.
- Near-axis cells lose the C-band geometries first; the axial vertical
  probe of the constrained variant rests on the two NISAR looks alone
  and is not a joint-sensor result.
- Closure is strongly position-dependent: it is smallest in the far
  field and grows in the deforming area, where window differences and
  steep gradients dominate; see the regional closure line above.
- Formal sigmas propagate the far-field dataset scatter only.
- The dike frame depends on the chosen strike; the along-strike
  component is the weakly observed one, and a few degrees of strike
  moves a few percent of U⊥ into U∥ (see the sensitivity line).
- The decay normalization assumes the post-event transient shares the
  co-event dike pattern; it is a model choice (see LIMITATIONS.md).
