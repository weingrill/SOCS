M6: The Butterfly Cluster
UPDATE clusters SET notes = 'M 6' WHERE name = 'NGC 6405';
UPDATE clusters SET notes = 'M 7' WHERE name = 'NGC 6475';
UPDATE clusters SET notes = 'M 11' WHERE name = 'NGC 6705';
UPDATE clusters SET notes = 'M 18' WHERE name = 'NGC 6613';
UPDATE clusters SET notes = 'M 21' WHERE name = 'NGC 6531';
UPDATE clusters SET notes = 'M 23' WHERE name = 'NGC 6494';
UPDATE clusters SET notes = 'M 25' WHERE name = 'IC 4725 ';
UPDATE clusters SET notes = 'M 26' WHERE name = 'NGC 6694';
UPDATE clusters SET notes = 'M 29' WHERE name = 'NGC 6913';
UPDATE clusters SET notes = 'M 34' WHERE name = 'NGC 1039';
UPDATE clusters SET notes = 'M 35' WHERE name = 'NGC 2168';
UPDATE clusters SET notes = 'M 36' WHERE name = 'NGC 1960';
UPDATE clusters SET notes = 'M 37' WHERE name = 'NGC 2099';
UPDATE clusters SET notes = 'M 38' WHERE name = 'NGC 1912';
UPDATE clusters SET notes = 'M 39' WHERE name = 'NGC 7092';
UPDATE clusters SET notes = 'M 41' WHERE name = 'NGC 2287';
UPDATE clusters SET notes = 'M 44' WHERE name = 'NGC 2632';
UPDATE clusters SET notes = 'M 45' WHERE name = 'Melotte 22';
UPDATE clusters SET notes = 'M 46' WHERE name = 'NGC 2437';
UPDATE clusters SET notes = 'M 47' WHERE name = 'NGC 2422';
UPDATE clusters SET notes = 'M 48' WHERE name = 'NGC 2548';
UPDATE clusters SET notes = 'M 50' WHERE name = 'NGC 2323';
UPDATE clusters SET notes = 'M 52' WHERE name = 'NGC 7654';
UPDATE clusters SET notes = 'M 67' WHERE name = 'NGC 2682';
UPDATE clusters SET notes = 'M 93' WHERE name = 'NGC 2447';
UPDATE clusters SET notes = 'M 103' WHERE name = 'NGC 581';

UPDATE clusters SET observed= TRUE WHERE name in ('NGC 2422', 'NGC 2323', 'NGC 2548', 'NGC 6633', 'NGC 6709','NGC 2281','NGC 2682', 'NGC 1528', 'NGC 6940');