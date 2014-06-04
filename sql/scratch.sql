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

