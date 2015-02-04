CREATE TABLE m48stars (
    starid character varying(25) NOT NULL,
    bv real,
    vmag real DEFAULT 0,
    vmag_err real,
    bmag real DEFAULT 0,
    bmag_err real,
    period real,
    period_err real,
    theta real,
    amp real,
    amp_err real,
    nv integer DEFAULT 0,
    nb integer DEFAULT 0,
    ra double precision,
    "dec" double precision,
    coord point,
    simbad character varying(25),
    good boolean,
    freq real,
    s1 real,
    c1 real,
    s2 real,
    c2 real,
    member boolean,
    clean_period real,
    clean_amp real,
    clean_sigma real,
    pman real,
    quality character varying(5),
    notes character varying(32),
    tab integer
);


select object, filter, expt
from frames
where object like 'M 48 BVI %'
group by object, filter, expt
order by object, filter, expt;

select object, min(datesend), max(datesend)
from frames where object like 'M 48 BVI %'
group by object;

select object, filter, count(objid) from frames where object like 'M 48 BVI %' group by object, filter order by object, filter;

SELECT phot.objid, mag_auto-corr 
FROM phot, frames
WHERE object like 'M 48 BVI%'
AND phot.objid=frames.objid
AND filter='V'
AND flags<8
AND point(123.326366,-5.808771) <@ circle(phot.coord,1./3600.)
ORDER BY abs(corr)
LIMIT 5;

SELECT count(objid)
FROM frames
WHERE object like 'M 48 BVI%%'
AND frames.good=True
AND filter='V';

SELECT objid,object, filter, good from frames
WHERE object like 'M 48 BVI%%';

