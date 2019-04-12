SELECT source_id, ra, dec, phot_g_mean_mag AS gmag, pmra, pmdec, 1000./parallax AS dist
FROM gaiadr2
WHERE CIRCLE(POINT(82.16, 35.8483), 1.0) @> CIRCLE(POINT(ra, dec), 0.0) AND phot_g_mean_mag < 18.5 AND phot_bp_mean_mag-phot_rp_mean_mag > 0.0 AND phot_bp_mean_mag-phot_rp_mean_mag < 2.5 AND 1000./parallax > 700 AND 1000./parallax < 2100;