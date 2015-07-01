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

UPDATE m48stars set good=NULL, provisional=NULL;
update m48stars set good=true where tab in (284, 287, 303, 313, 332, 336, 343, 363, 368, 386, 407, 418, 421, 425, 437, 467, 482, 485, 501, 507, 511, 517, 536, 540, 555, 556, 602, 620, 633, 649, 652, 657, 674, 679, 699, 713, 796, 807, 862, 864, 872, 898, 909, 920, 921, 923, 931, 935, 937, 969, 974, 975, 1227, 1455, 1456, 1605, 1711, 1744, 2010, 2071, 2285, 2346, 2632);
update m48stars set provisional=true where tab in (316, 562, 752, 772, 954, 1096, 1162);

--- update m48stars set provisional=false where good;
update m48stars set good=true where provisional or not provisional;
select count(tab),good,provisional from m48stars group by good, provisional;

