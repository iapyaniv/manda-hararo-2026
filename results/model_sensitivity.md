# Uniform-dike model sensitivity

These runs change one processing choice at a time. They are **not**
confidence intervals and do not account for spatially correlated atmosphere,
faulting above the dike, nonuniform opening, or elastic-structure error.

| variant | n | RMS (cm) | strike | length (km) | depth (km) | opening (m) | volume (km³) |
|---|---:|---:|---:|---:|---:|---:|---:|
| default | 7717 | 3.05 | 152.6° | 17.15 | 1.49–4.99 | 2.81 | 0.169 |
| 0.96 km blocks | 13575 | 3.06 | 152.2° | 17.19 | 1.53–4.87 | 2.92 | 0.168 |
| 1.92 km blocks | 3357 | 2.96 | 152.2° | 17.51 | 1.30–5.50 | 2.36 | 0.174 |
| 25 km far field | 7717 | 3.07 | 152.6° | 17.14 | 1.51–4.93 | 2.85 | 0.167 |
| 32 km far field | 7717 | 3.03 | 152.6° | 17.17 | 1.47–5.07 | 2.75 | 0.170 |
| LOS height 0 m | 7717 | 3.05 | 152.6° | 17.15 | 1.49–4.99 | 2.81 | 0.169 |
| LOS height 1000 m | 7717 | 3.05 | 152.6° | 17.15 | 1.50–4.99 | 2.81 | 0.168 |
| equal geometry weight | 7717 | 3.07 | 152.3° | 17.43 | 1.45–4.96 | 2.80 | 0.171 |
| fit within 30 km | 2919 | 4.28 | 152.7° | 17.33 | 1.30–6.00 | 2.31 | 0.188 |

Envelope across these runs:

- length 17.1–17.5 km
- top depth 1.3–1.5 km; bottom depth 4.9–6.0 km
- opening 2.31–2.92 m
- geometric volume 0.167–0.188 km³

The volume varies less than depth and opening separately. A statistical uncertainty
analysis requires a spatial covariance model and comparison with nonuniform-opening
and faulting models.
