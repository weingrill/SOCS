CREATE TABLE ngc6633 (
    starid character varying(25) NOT NULL,
    bv real,
    vmag real DEFAULT 0,
    vmag_err real,
    nv integer DEFAULT 0,
    bmag real DEFAULT 0,
    bmag_err real,
    nb integer DEFAULT 0,
    period real,
    period_err real,
    theta real,
    amp real,
    ra double precision,
    "dec" double precision,
    coord point,
    ref character varying(10),
    good boolean,
    freq real,
    member boolean,
    clean_period real,
    clean_amp real,
    clean_sigma real,
    notes character varying(32),
    tab integer
);

CREATE TABLE ngc6633ref (
    starid character varying(10) NOT NULL,
    ra double precision,
    "dec" double precision,
    dx real,
    dy real,
    x real,
    y real,
    bmag real,
    bsigma real,
    nb integer,
    vmag real,
    vsigma real,
    nv integer,
    imag real,
    isigma real,
    ni integer,
    coord point
);

COMMENT ON TABLE ngc6633ref IS 'Stetson 2000PASP..112..925S';
CREATE INDEX idx_ngc6633ref_coord ON ngc6633ref USING GIST (circle(coord,0));

INSERT INTO ngc6633ref (starid, ra, dec, dx, dy, x, y, coord) 
VALUES ('Reference', 276.81317437092,  +06.53027772692, 0.00, 0.00, 5170.5, 4386.5, point(276.81317437092,  +06.53027772692)):

update ngc6633ref set bmag=NULL WHERE bmag>99;
