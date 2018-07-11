import astropy
import jplephem
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric
t = Time("2014-09-22 23:22")
with solar_system_ephemeris.set('de432s'):
    jup = get_body_barycentric(['saturn','moon'], t)

