SELECT id from matched
FROM matched, frames, phot
WHERE frames.objid = matched.objid
 AND phot.objid = matched.objid
 AND phot.star = matched.star
 AND frames.object LIKE 'NGC 1647 field 1'

CREATE table ngc1647stars AS
SELECT DISTINCT matched.id
FROM matched, frames, phot
WHERE frames.object LIKE 'NGC 1647 field 1'
 AND frames.objid = phot.objid
 AND phot.objid = matched.objid
 AND phot.star = matched.star;

select objid, round(100.0*matched/stars,3) "found", stars,expt
from frames 
WHERE frames.object LIKE 'NGC 1647 field 1'
order by found desc
limit 20;

CREATE VIEW ngc1647f1 AS
SELECT objid
FROM frames
WHERE object like 'NGC 1647 field 1' 
 AND frames.expt >= 60 AND filter like 'rp' 
ORDER by objid;

update frames
set good=False
where 
object like 'NGC 1647 field 1'
and expt=62.5
and good is NULL;

select * from phot where objid='20101022A-0000-0003';

SELECT id, mag_isocor 
FROM matched, phot
WHERE phot.objid='20101022A-0000-0003'
 AND phot.objid = matched.objid
 AND phot.star = matched.star;

SELECT ngc1647ref.id, ngc1647ref.mag "ref", ngc1647ref.std, phot.mag_isocor "mag"
FROM matched, phot, ngc1647ref
WHERE phot.objid='20101011A-0004-0012'
 AND phot.objid = matched.objid
 AND phot.star = matched.star
 AND ngc1647ref.id=matched.id;

drop view ngc1647ref;
create table ngc1647ref as 
SELECT id, avg(mag_isocor) "mag", stddev_samp(mag_isocor) "std"
FROM matched, phot
WHERE phot.objid like '20101022A-0000-000%'
 AND phot.objid = matched.objid
 AND phot.star = matched.star
group by id;
grant select on ngc1647ref to sro;
grant select on ngc1647ref to jwe;

select * from matched where objid like '20101031A-0012-0002';

select matched.id, phot.flags
from matched, phot 
where matched.objid like '20101031A-0012-0002'
 and matched.objid=phot.objid
 and matched.star=phot.star;

CREATE TABLE ngc1647stars (
 id varchar(25),
 vmag real,
 bv real,
 period real,
 period_err real,
 PRIMARY KEY (id)
);

select object, expt, count(objid) 
from frames 
where object like '%1647%' 
 and filter like 'rp'
group by object, expt;

update frames
set good=true
where

insert into ngc1647stars (id, ra, dec)
select distinct matched.id "id", ppmxl.ra "ra", ppmxl.dec "dec"
from matched, phot, ppmxl
where matched.objid like '20101031A-0012-0002'
 and (matched.objid,matched.star)=(phot.objid,phot.star)
 and (matched.id = ppmxl.id);

insert into ngc1647stars
select distinct matched.id 
from matched, phot, ppmxl
where matched.objid like '20101031A-0012-0002'
 and (matched.objid,matched.star)=(phot.objid,phot.star)
 and (matched.id = ppmxl.id);

SELECT distinct matched.id, ppmxl.ra, ppmxl.dec
FROM frames, phot, matched, ppmxl
WHERE frames.object like 'NGC 1647 bv %'
AND frames.filter like 'V'
AND frames.objid = phot.objid
AND (matched.objid,matched.star) = (phot.objid,phot.star)
AND matched.id = ppmxl.id;

SELECT distinct matched.id, ppmxl.ra, ppmxl.dec
FROM frames, phot, matched, ppmxl
WHERE frames.object like 'NGC 1647 field %'
AND frames.filter like 'rp'
AND frames.objid = phot.objid
AND (matched.objid,matched.star) = (phot.objid,phot.star)
AND matched.id = ppmxl.id;

UPDATE ngc1647stars
SET (ra,dec) = (ppmxl.ra,ppmxl.dec)
FROM ppmxl
WHERE ngc1647stars.id = ppmxl.id

