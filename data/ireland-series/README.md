# Digitized site series from B. Ireland's Sentinel-1 time-series figure

`ireland_079D_sites_digitized.csv` holds the six Sentinel-1 descending
(079D_07694_131313) LOS series of Ben Ireland's (Bristol) time-series figure
of 2026, read off its raster panels: the
marker centres were calibrated from the axis ticks and assigned to the public
079D acquisition dates within 20 days (295 to 306 of 313 dates per site;
one pixel is 0.6 to 1.6 mm, the median date offset 2 to 4 days). Values are
LOS in millimetres, positive toward the satellite, relative to the first date
(23 October 2014) in the original processing reference. The "Graben
subsidence" series is the graben site minus the flank site 8.5 km to the
north-east, as in the figure.

This is a draft stand-in for `mh26.figure1` until the series themselves are
supplied; it is not a measurement to be cited.

Columns: site, lon, lat, epoch, los_mm, x_px, y_px, dt_days, blob_markers.
