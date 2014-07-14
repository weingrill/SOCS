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
