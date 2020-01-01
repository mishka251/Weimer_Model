import numpy as np
from w05sc import Calculator

nlat = 97
nlon = 25

by = 0.
bz = -5.
tilt = 0.

swvel = 450.
swden = 9.
file_path = "./"

fill = 1.e36

epot = np.zeros((nlon, nlat), np.float)
fac = np.zeros((nlon, nlat), np.float)

calc = Calculator()
calc.setmodel(by, bz, tilt, swvel, swden, file_path, 'epot')

mlat = [-90., -88.1238292398491, -86.2386359278657, -84.3344382773342,
        -82.4013318763435, -80.4295344892688, -78.4094552099168,
        -76.331796630125, -74.1876988925388, -71.9689341802758,
        -69.6681589022773, -67.2792279882741, -64.7975706790533,
        -62.2206194320588, -59.5482728298363, -56.7833601290164,
        -53.9320608459732, -51.0042204168578, -48.0134966005524,
        -44.9772754602266, -41.916313892128, -38.8540980954293,
        -35.8159497801506, -32.8279553674349, -29.9158266703621,
        -27.1038148776609, -24.4137889090065, -21.8645574169981,
        -19.4714697638694, -17.2462861630082, -15.1972697734841,
        -13.3294282264571, -11.6448185129562, -10.142824406667,
        -8.82031765103987, -7.67162666281269, -6.68827297583048,
        -5.85851734698832, -5.16689314460211, -4.5940469432968,
        -4.11722526306697, -3.71151170575937, -3.35148255039153,
        -3.01257883277328, -2.67136426606314, -2.3036287214954,
        -1.87754943767857, -1.32687203939232, -7.72840966450717e-08,
        1.32687203939232, 1.87754943767857, 2.3036287214954, 2.67136426606314,
        3.01257883277328, 3.35148255039153, 3.71151170575936, 4.11722526306697,
        4.59404694329679, 5.16689314460211, 5.85851734698832, 6.68827297583048,
        7.67162666281268, 8.82031765103987, 10.142824406667, 11.6448185129562,
        13.3294282264571, 15.1972697734841, 17.2462861630082, 19.4714697638694,
        21.8645574169981, 24.4137889090064, 27.1038148776609, 29.9158266703621,
        32.8279553674348, 35.8159497801506, 38.8540980954293, 41.916313892128,
        44.9772754602266, 48.0134966005524, 51.0042204168578, 53.9320608459731,
        56.7833601290163, 59.5482728298363, 62.2206194320588, 64.7975706790533,
        67.2792279882741, 69.6681589022773, 71.9689341802758, 74.1876988925387,
        76.331796630125, 78.4094552099168, 80.4295344892687, 82.4013318763434,
        84.3344382773342, 86.2386359278657, 88.123829239849, 90.]

mlt = [i for i in range(nlon)]
lon = [15 * i for i in range(nlon)]

for i in range(nlon):
    for j in range(nlon):
        epot[i][j] = calc.epotval(mlat[j], mlt[i], fill)

calc.setmodel(by, bz, tilt, swvel, swden, file_path, 'bpot')

for i in range(nlon):
    for j in range(nlon):
        fac[i][j] = calc.mpfac(mlat[j], mlt[i], fill)

datfile = 'wei05sc_epot_f90.dat '
file = open(file=datfile, mode="w")
file.write(str(nlon) + "  " + str(nlat))
np.savetxt(file, mlt)  # file.write(str(mlt))
np.savetxt(file, mlat)  # file.write(str(mlat))
np.savetxt(file, epot)  # file.write(str(epot))
file.close()
print(f"('Wrote ascii file '{datfile})")

datfile = 'wei05sc_fac_f90.dat '
file = open(datfile, mode="w")
file.write(f"{nlon}, {nlat}")
np.savetxt(file, mlt)  # file.write(str(mlt))
np.savetxt(file, mlat)  # file.write(str(mlat))
np.savetxt(file, fac)  # file.write(str(fac))
file.close()
print(f"('Wrote ascii file '{datfile})")
