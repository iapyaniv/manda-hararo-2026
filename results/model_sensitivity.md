# Uniform-dike model sensitivity

These runs change one processing choice at a time. They are **not**
confidence intervals and do not account for spatially correlated atmosphere,
faulting above the dike, nonuniform opening, or elastic-structure error.

| variant | n | RMS (cm) | strike | length (km) | depth (km) | opening (m) | volume (km³) |
|---|---:|---:|---:|---:|---:|---:|---:|
| default | 7717 | 2.91 | 153.0° | 17.31 | 1.30–5.94 | 2.32 | 0.187 |
| 0.96 km blocks | 13575 | 2.92 | 152.6° | 17.35 | 1.33–5.78 | 2.40 | 0.186 |
| 1.92 km blocks | 3357 | 2.83 | 152.5° | 17.71 | 1.14–6.46 | 2.05 | 0.193 |
| 25 km far field | 7717 | 2.91 | 153.0° | 17.32 | 1.30–5.95 | 2.32 | 0.187 |
| 32 km far field | 7717 | 2.91 | 153.0° | 17.31 | 1.30–5.94 | 2.33 | 0.187 |
| LOS height 0 m | 7717 | 2.91 | 153.0° | 17.31 | 1.30–5.94 | 2.32 | 0.187 |
| LOS height 1000 m | 7717 | 2.91 | 153.0° | 17.32 | 1.30–5.94 | 2.33 | 0.187 |
| equal geometry weight | 7717 | 2.93 | 152.6° | 17.57 | 1.30–5.77 | 2.38 | 0.187 |
| fit within 30 km | 2919 | 4.21 | 153.0° | 17.54 | 1.18–6.86 | 2.08 | 0.208 |

Envelope across these runs:

- length 17.3–17.7 km
- top depth 1.1–1.3 km; bottom depth 5.8–6.9 km
- opening 2.05–2.40 m
- geometric volume 0.186–0.208 km³

The volume varies less than depth and opening separately. A statistical uncertainty
analysis requires a spatial covariance model and comparison with nonuniform-opening
and faulting models.
