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

