# Dike-induced stress change on the saved fault surfaces

- source fit: `results/distributed_fit_mapped.json` (SHA-256 `80a46f46dff248c6b173d4fb1bff13fc4baaff29dc505336b0d752cfa0ccbfea`)
- source: distributed opening with trace-column geometry; fault slip is excluded from the stress source
- receivers: patch centres from the fit's saved trace columns; dip-specific truncation at the dike plane; 400 m dike standoff
- spatial output: normal and Coulomb stress on a 250-m-depth by saved-column grid for the 65° receiver case

## Median normal-stress change by receiver depth (dip 65°; unclamping positive)

- east: 0.7 km: +5.2 MPa; 2.0 km: -2.2 MPa; 3.3 km: -10.1 MPa (overall +0.4 MPa, 50% of receivers unclamped)

## Receiver-dip sensitivity

- dip 55°: east +1.4 / +1.1 MPa
- dip 65°: east +0.4 / -0.1 MPa
- dip 75°: east -0.7 / -0.7 MPa

Normal stress varies with receiver depth. The shallowest saved row is mainly unclamped, while most deeper saved rows are clamped in the 65-degree case. This result is conditional on the saved receiver surfaces. The distributed opening starts at the surface and has no subsurface hard top, so this report makes no conclusion relative to the uniform rectangle's fitted top.

## Scope

- Static elastic half-space stress change from the inverted opening
  only; no fault interaction, pore pressure, or topography.
- Source patches and receiver trace columns come from the saved fit
  geometry. Current source-code geometry constants are not used.
- Receiver dip changes are sensitivity cases. The source geometry
  does not change between these cases.
