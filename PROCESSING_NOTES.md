# Processing notes

This file records processing details that affect the numerical results.

## GUNW extraction

The GUNW HDF5 coherence dataset is named `coherenceMagnitude`. The extraction
code maps source dataset names to the shorter NPZ keys and checks the required
datasets, raster shapes, projection, metadata, and output file before use.

## Valid-sample mask

Coherence and connected-component screens do not remove every sample that the
GUNW product marks as invalid. All analysis paths and GeoTIFF exports therefore
apply the reference and secondary valid-sample flags before further screening.

Without this mask, the main ascending fit pair contains 21,183 invalid pixels
(2.8% of the earlier analysis mask). Each track-100 pair contains about 45,000
such pixels (7.3–7.6%). The main descending fit pair contains none. These
pixels produce a diagonal subswath-edge feature in the ascending field.

No block is rejected on the basis of model residuals.

## Inversion scaling

The inversion contains kilometre-scale position and length parameters and
metre-scale opening and offset parameters. Optimization uses dimensionless
scaled variables and explicit finite-difference steps. All six starting points
converge to the same solution basin. The best optimality is 0.0043, and the
maximum length difference among the six solutions is less than 1 m.

## Sentinel-1 coherence

The local LiCSAR coherence raster uses UInt8 values from 0 to 255. The
comparison normalizes this raster to 0–1 before it applies the 0.25 threshold.
It also requires identical grids for all five input rasters and checks the LOS
vector norms.

## Parameter sensitivity

The formal covariance in `results/okada_fit.json` assumes independent white
noise and is only a numerical diagnostic. Nine processing-choice runs vary
block size, far-field radius, LOS interpolation height, geometry weighting,
and fit radius. Their range describes processing sensitivity, not statistical
uncertainty.

The valid-sample and optimization changes alter the fitted geometry:

| parameter | earlier processing | current processing |
|---|---:|---:|
| length | 18.1 km | 17.15 km |
| depth | 0.8–7.0 km | 1.49–4.99 km |
| uniform opening | 1.58 m | 2.81 m |
| geometric volume | 0.18 km³ | 0.169 km³ |
| NISAR RMS | 3.8 cm | 3.05 cm |

All files in `results/` use the current processing.

## Archive verification

The core, pre-event, and complete data sets have separate SHA-256 manifests.
The core archive can therefore be verified without the optional pre-event NPZ
files. Release archives use the source commit time, record the source commit
hash, and pass a ZIP CRC check before replacement of an existing archive.

## Additional analysis

The following work is not part of the current results:

- estimation with a spatial covariance model;
- comparison with alternative source geometries;
- comparison of the post-event pair with additional viewing geometries;
- integration with long-term deformation and local seismic data.

The validator checks data schema, coordinates, masks, LOS vectors, raster
alignment, numerical ranges, source manifests, and SHA-256 hashes. The tests
cover the validity mask, coherence encoding, block averaging, and LOS-vector
interpolation.
