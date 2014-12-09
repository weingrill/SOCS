INSERT INTO ngc2281stars (id, ra, dec)
SELECT DISTINCT matched.id, ppmxl.ra, ppmxl.dec
FROM frames, phot, matched, ppmxl
WHERE frames.object LIKE 'NGC 2281 BVI %'
AND frames.filter LIKE 'V'
AND frames.objid = phot.objid
AND (matched.objid,matched.star) = (phot.objid,phot.star)
AND matched.id = ppmxl.id;

SELECT frames.hjd, phot.mag_auto, phot.magerr_auto, phot.flags
FROM frames, matched, phot
WHERE id LIKE 'PPMXL3365775734594919249'
AND filter LIKE 'V'
AND frames.objid = matched.objid
AND (phot.objid,phot.star) = (matched.objid,matched.star)
ORDER BY hjd;

CREATE TABLE ngc2281ref (
 id varchar(8),
 ra double precision,
 dec double precision,
 s char,
 ucac2 int,
 tyc2 varchar(11),
 vmag real,
 remarks varchar(50),
 bv real,
 coord point,
 PRIMARY KEY (id)
);

alter table ngc2281ref add column bmag real;

update ngc2281ref set bmag=bv+vmag;
update ngc2281ref set coord = point(ra,dec);

SELECT phot.mag_auto, phot.magerr_auto, vmag, bv
 FROM frames, phot, ngc2281ref
 WHERE frames.objid like '20140305A-0098-0005'
 AND frames.objid=phot.objid
 AND frames.filter='V'
 AND circle(ngc2281ref.coord,3./3600.) @> circle(phot.coord,.0)
 AND bv>0.0;

SELECT frames.hjd, phot.mag_auto, phot.magerr_auto
FROM frames, matched, phot
WHERE matched.id LIKE 'PPMXL2860071617461126379'
AND frames.object like 'M 48 rot%%'
AND filter LIKE 'V'
AND frames.objid = matched.objid
AND (phot.objid,phot.star) = (matched.objid,matched.star)
AND phot.flags<8
ORDER BY hjd;


SELECT id 
FROM matched
WHERE (matched.objid, matched.star) = ('20140303A-0074-0013',2052);

SELECT frames.hjd, phot.mag_auto-corr, phot.magerr_auto, phot.flags
FROM frames, matched, phot
WHERE matched.id LIKE 'PPMXL2859922316456747397'
AND frames.object like 'M 48 rot%%'
AND filter LIKE 'V'
AND frames.good
AND NOT corr IS NULL
AND frames.objid = matched.objid
AND (phot.objid,phot.star) = (matched.objid,matched.star)
AND phot.flags<4
ORDER BY hjd;

COMMENT ON COLUMN clusters.name IS 'Cluster name';
COMMENT ON COLUMN clusters.ra IS '[deg] Right ascension (J2000.0)';
COMMENT ON COLUMN clusters.dec IS '[deg] Declination (J2000.0)';
COMMENT ON COLUMN clusters.class IS '[*] Flag for classification of the cluster';
COMMENT ON COLUMN clusters.diam IS '[arcmin] Apparent diameter in arcmin';
COMMENT ON COLUMN clusters.d IS '[pc] Distance';
COMMENT ON COLUMN clusters.ebv IS 'Colour excess in BV';
COMMENT ON COLUMN clusters.logage IS '[yr] Age (in log t)';
COMMENT ON COLUMN clusters.pmra IS '[mas/yr] Mean proper motion of the cluster in mu_alpha.cos(delta), ICRS';
COMMENT ON COLUMN clusters.epmra IS '[mas/yr] Standard deviation in pmRA, ICRS';
COMMENT ON COLUMN clusters.pmdec IS '[mas/yr] Mean proper motion of the cluster in mu_delta, ICRS';
COMMENT ON COLUMN clusters.epmdec IS '[mas/yr] Standard deviation in pmRA and pmDE, ICRS';
COMMENT ON COLUMN clusters.nc IS 'Estimated number of members in the cluster';
COMMENT ON COLUMN clusters.ref1 IS 'Source of the mean proper motion determination';
COMMENT ON COLUMN clusters.rv IS '[km/s] Radial Velocity ';
COMMENT ON COLUMN clusters.erv IS '[km/s] Error in Radial Velocity';
COMMENT ON COLUMN clusters.n IS 'Number of stars used to determine Radial Velocity';
COMMENT ON COLUMN clusters.ref2 IS 'Source of the mean radial velocity determination';
COMMENT ON COLUMN clusters.me IS 'Metallicity';
COMMENT ON COLUMN clusters.eme IS 'Error in Metallicity';
COMMENT ON COLUMN clusters.nme IS 'Number of stars used to determine Metallicity';
COMMENT ON COLUMN clusters.trtyp IS 'Trumpler Type determined in the DSS inspection';

