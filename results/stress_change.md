# Dike-induced stress change on the saved fault surfaces

- source fit: `results/distributed_fit_mapped.json` (SHA-256 `8c287924a3eac112f4073cda177240fc86c6b45d46bb9088210a1d0bc8308f32`)
- source: distributed opening with trace-column geometry; fault slip is excluded from the stress source
- receivers: patch centres from the fit's saved trace columns; dip-specific truncation at the dike plane; 400 m dike standoff
- spatial output: normal and Coulomb stress on a 250-m-depth by saved-column grid for the 65° receiver case

## Median normal-stress change by receiver depth (dip 65°; unclamping positive)

- east: 0.7 km: +6.1 MPa; 2.0 km: -3.1 MPa; 3.3 km: -8.3 MPa (overall +1.0 MPa, 52% of receivers unclamped)

## Receiver-dip sensitivity

- dip 55°: east +1.6 / +1.0 MPa
- dip 65°: east +1.0 / +0.2 MPa
- dip 75°: east -0.3 / -0.3 MPa

Normal stress varies with receiver depth. The shallowest saved row is mainly unclamped, while most deeper saved rows are clamped in the 65-degree case. This result is conditional on the saved receiver surfaces. The distributed opening starts at the surface and has no subsurface hard top, so this report makes no conclusion relative to the uniform rectangle's fitted top.

## Scope

- Static elastic half-space stress change from the inverted opening
  only; no fault interaction, pore pressure, or topography.
- Source patches and receiver trace columns come from the saved fit
  geometry. Current source-code geometry constants are not used.
- Receiver dip changes are sensitivity cases. The source geometry
  does not change between these cases.
