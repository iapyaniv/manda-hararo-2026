# Analysis guide

This repository contains NISAR displacement fields, a uniform-opening source
model, a Sentinel-1 comparison, a pre-event probe history, a post-event
pattern projection, and exact source-product identifiers for the July 2026
Manda Hararo intrusion.

## Package contents

- `results/los_panels.png` and `results/decomposition.png` show the NISAR LOS
  fields and the two-geometry quasi-East/quasi-Up decomposition.
- `results/okada_fit_panels.png` and `results/okada_fit.json` contain the
  uniform-opening dike fit.
- `results/model_sensitivity.md` reports the effect of nine processing choices.
- `results/s1_crosscheck_panels.png` compares the NISAR-derived model with a
  Sentinel-1 LiCSAR interferogram.
- `results/probe_history.md` reports 16 pre-event windows and the co-event and
  post-event windows.
- `results/continuation_panels.png` compares the 19 July–12 August descending
  pair with the co-event dike pattern.
- `granules/` records the NISAR and LiCSAR source products.
- `data/` contains archive manifests and instructions for the data files.

The core data archive supports all reported results. The pre-event archive is
needed only to regenerate the 16-window history from the NPZ layers. The core
archive already contains the corresponding GeoTIFFs.

## Current source model

The reference model is a vertical, uniform-opening rectangle in a homogeneous
elastic half-space. Its fitted parameters are:

- strike: 152.6°;
- length: 17.15 km;
- depth: 1.49–4.99 km;
- opening: 2.81 m;
- geometric volume: 0.169 km³;
- NISAR RMS residual: 3.05 cm.

The processing-choice runs span 0.167–0.188 km³. This range is not a
confidence interval. The separate depth and opening estimates are more
sensitive to processing choices than their product.

The model uses the GUNW valid-sample mask and scaled optimization variables.
The effects of these choices are described in `PROCESSING_NOTES.md`.

## Post-event pair

For the descending 19 July–12 August pair, the projection onto the co-event
dike pattern is +0.107. The processing-choice range is +0.100 to +0.123. The
field also contains a broad residual and a fitted +3.3 cm subset offset. One
descending geometry cannot distinguish renewed opening from atmospheric or
other deformation signals. The coefficient is therefore not converted to an
added-volume estimate.

## Scope

- Long-term regional InSAR time series are not included.
- Local seismic waveform and catalogue data are not included.
- The source model does not include fault slip, topography, elastic layering,
  or a coupled reservoir.
- No seismic-to-geodetic moment ratio is calculated.
- No relation to other Afar events is tested.
- Public reuse terms have not been selected.

See `LIMITATIONS.md` before interpreting the model parameters or the
post-event projection.

## Suggested processing order

```bash
python -m unittest discover -s tests -v
python -m mh26.validate --profile core
python -m mh26.los_fields
python -m mh26.invert
python -m mh26.sensitivity
python -m mh26.s1_check
python -m mh26.history
python -m mh26.continuation
```

The Markdown reports and JSON files in `results/` contain the same numerical
outputs in human-readable and machine-readable forms.
