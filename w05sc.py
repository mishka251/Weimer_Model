from math import sqrt, atan2, radians, degrees, sin, exp, cos, atan, log, asin, factorial
from typing import List, Tuple

from reader import Reader
import numpy as np

reader: Reader = Reader()

bndyfitr: float
esphc: List[float] = [0] * reader.csize
bsphc: List[float] = [0] * reader.csize
tmat: np.ndarray = np.zeros(3, 3)
ttmat: np.ndarray = np.zeros(3, 3)
#
mxtablesize: int = 200

#
plmtable: np.ndarray = np.zeros(mxtablesize, reader.csize)
colattable: np.ndarray = np.zeros(mxtablesize)
#
nlms: np.ndarray = np.zeros(reader.csize)


# def checkinputs(lat, mlt, inside, phir, colat):
#     pass


def lngamma(xx: float) -> float:
    # This is an f90-python translation from C code copied from
    # www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
    x: float
    y: float
    tmp: float
    ser: float

    cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
           -0.5395239384953e-5]

    y = xx
    x = xx
    tmp = x + 5.5
    tmp = tmp - (x + 0.5) * log(tmp)
    ser = 1.000000000190015

    for j in range(len(cof)):
        y = y + 1
        ser = ser + cof[j] / y

    return -tmp + log(2.5066282746310005 * ser / x)


def value_locate(vec: List[float], val: float) -> int:
    if val < vec[0]:
        return -1
    if val > vec[-1]:
        return len(vec)
    for i in range(len(vec) - 1):
        if vec[i] <= val <= vec[i + 1]:
            return i


def interpol_quad(v: List[float], x: List[float], u: List[float]) -> List[float]:
    nv: int = len(v)
    nx: int = len(x)
    nu: int = len(u)

    assert nx == nv, f"('>>> interpol_quad: nx /= nv: nx='{nx}' nv='{nv})"
    p: List[float] = []
    for i in range(nu):
        ix: int = value_locate(x, u[i])
        if ix <= 1 or ix >= nx:  # ! bug fix by btf 12/23/09
            p.append(0.)
            continue  # cycle                       ! bug fix by btf 12/23/09
        # endif
        x1 = x[ix]
        x0 = x[ix - 1]
        x2 = x[ix + 1]
        p_i = v[ix - 1] * (u[i] - x1) * (u[i] - x2) / ((x0 - x1) * (x0 - x2)) + \
              v[ix] * (u[i] - x0) * (u[i] - x2) / ((x1 - x0) * (x1 - x2)) + \
              v[ix + 1] * (u[i] - x0) * (u[i] - x1) / ((x2 - x0) * (x2 - x1))

        p.append(p_i)
    return p


def dorotation(latin: float, lonin: float) -> Tuple[float, float]:
    """

    :param latin:
    :param lonin:
    :return:latout,lonout
    """
    # latr: float
    # lonr: float
    # stc: float
    # ctc: float
    # sf: float
    # cf: float
    # a: float
    # b: float
    pos = np.zeros(3)

    latr: float = radians(latin)
    lonr: float = degrees(lonin)
    stc: float = sin(latr)
    ctc: float = cos(latr)
    sf: float = sin(lonr)
    cf: float = cos(lonr)
    a: float = ctc * cf
    b: float = ctc * sf

    for i in range(3):
        pos[i] = tmat[0, i] * a + tmat[1, i] * b + tmat[2, i] * stc

    latout: float = degrees(asin(pos[0]))
    lonout: float = degrees(atan2(pos[1], pos[0]))

    return latout, lonout


def checkinputs(lat: float, mlt: float) -> Tuple[int, float, float]:
    """

    :param lat:
    :param mlt:
    :return:  inside: int, phir: float, colat: float
    """
    lon: float
    tlat: float
    tlon: float
    radii: float

    lon: float = mlt * 15.
    tlat, tlon = dorotation(lat, lon)
    radii = 90. - tlat
    inside: int = 1 if (radii <= bndyfitr) else 0
    phir: float = radians(tlon)
    colat: float = radii

    return inside, phir, colat


def nkmlookup(k: int, m: int, th0: float) -> float:
    if th0 == 90.:
        return float(k)

    kk: int = k + 1
    mm: int = m + 1

    th0a: float = th0

    res: List[float] = []

    if kk >= reader.maxk_scha:
        print(f"('>>> nkmlookup: kk > maxk: kk='{kk}' maxk='{reader.maxk_scha}")
        res = interpol_quad(reader.allnkm[reader.maxk_scha - 1][mm], reader.th0s, [th0a])
    if mm >= reader.maxm_scha:
        print(f"('>>> nkmlookup: mm > maxm: kk='{kk}' maxm='{reader.maxm_scha}")
        res = interpol_quad(reader.allnkm[kk][reader.maxm_scha - 1], reader.th0s, [th0a])
    if th0 < reader.th0s[0]:
        print(f"('>>> nkmlookup: th0 < th0s(1): th0='{th0}' th0s(1)='{reader.th0s[0]}")

    res = interpol_quad(reader.allnkm[kk][mm], reader.th0s, [th0a])
    return res[0]


def km_n(m: int, rn: float) -> float:
    if m == 0:
        return 1

    return sqrt(2. * exp(lngamma(rn + m + 1.) - lngamma(rn - m + 1.))) / (2. ** m * factorial(m))


def pm_n(m: int, r: float, cth: List[float], tablesize: int) -> List[float]:
    a: List[float]
    if m == 0:
        a = [0]
    else:
        a = [sqrt(1 - cth[i] ** 2) ** m for i in range(1, tablesize + 1)]
    xn: float = r * (r + 1)

    x = [(1 - ct) / 2 for ct in cth]

    table: List[float] = a

    tmp: List[float] = [10000] * tablesize
    k = 0
    while max(tmp) > 1e-6:
        for i in range(tablesize):
            a[i] *= (x[i] * ((k + m - 1.) * (k + m) - xn) / (k * (k + m)))
            table[i] += a[i]
            pass

        k += 1

        for i in range(tablesize):
            div = abs(table[i])
            div = max(div, 1e-6)
            tmp[i] = abs(a[i]) / div
            pass
        pass

    ans = km_n(m, r)
    return [t * ans for t in table]


def scplm(index: int, colat: float, nlm: float) -> float:
    global plmtable, nlms, colattable
    # nlms:List[float] = [0]
    skip: bool = False
    th0 = bndyfitr
    if scplm.prevth0 != th0:
        scplm.tablesize = 3 * int(round((th0)))
        assert scplm.tablesize <= mxtablesize, f"('>>> tablesize > mxtablesize: tablesize={scplm.tablesize} mxtablesize={mxtablesize} tn0={th0}"

        colattable = [float(i - 1) * (th0 / float(scplm.tablesize - 1)) for i in range(scplm.tablesize)]
        cth = [cos(degrees(col)) for col in colattable]
        scplm.prevth0 = th0
        nlms = [0] * reader.csize

        for j in range(reader.csize):
            if skip:
                skip = False
                continue
                pass
            nlms.append(nkmlookup(reader.ls[j], reader.ms[j], th0))
            tmp: List[float] = pm_n(reader.ms[j], nlms[j], cth, scplm.tablesize)
            for _i in range(scplm.tablesize):
                plmtable[_i][index] = tmp[_i]
                pass
            skip = False

            if reader.ms[j] == 0 and reader.ab[j] > 0:
                plmtable[0][j + 1] = plmtable[0][j]
                nlms[j + 1] = nlms[j]
                skip = True
                pass
            pass  # end for
        pass  # endif

    nlm = nlms[index]
    colata = [colat]

    tmp = [plmtable[i][index] for i in range(scplm.tablesize)]

    out = interpol_quad(tmp, colattable[1: scplm.tablesize], colata)
    return out[1]


scplm.prevth0 = 1.e36
scplm.tablesize = 0


def mpfac(lat: float, mlt: float, fill: float) -> float:
    """

    :param lat:
    :param mlt:
    :param fill:
    :return: fac:float
    """
    ls = reader.ls
    ms = reader.ms
    ab = reader.ab

    m: int
    inside: int

    cfactor: float

    re: float
    z: float
    phir: float
    plm: float
    colat: float
    nlm: float = 0
    # pi:float

    re = 6371.2 + 110.  # km radius (allow default ht=110)

    inside, phir, colat = checkinputs(lat, mlt)

    if (inside == 0):
        return fill

    phim: List[float] = [phir, phir * 2]
    cospm: List[float] = [cos(phi) for phi in phim]
    sinpm: List[float] = [sin(phi) for phi in phim]

    z = 0.
    # jloop: do
    skip: bool = False
    for j in range(reader.csize):
        if skip:
            skip = False
            continue

        if ls[j] >= 11:
            break
        # m = ms[j]
        if ab[j] == 1:
            plm = scplm(j, colat, nlm)
            plm = plm * (nlm * (nlm + 1.))

            if ms[j] == 0:
                z = z - plm * bsphc[j]
            else:
                z = z - (plm * (bsphc[j] * cospm[ms[j]] + bsphc[j + 1] * sinpm[ms[j]]))
                skip = True

    pi = 4. * atan(1.)
    cfactor = -1.e5 / (4. * pi * re ** 2)  # convert to uA/m2
    z = z * cfactor
    return z


def epotval(lat: float, mlt: float, fill: float) -> float:
    """

    :param lat:
    :param mlt:
    :param fill:
    :param epot:
    :return: erot:float
    """
    # global inside, phir, colat
    # inside: int
    # j: int
    m: int
    mm: int

    z: float
    # phir: float
    plm: float
    # colat: float
    nlm: float

    # phim = np.zeros(2)
    # cospm = np.zeros(2)
    # sinpm = np.zeros(2)

    inside, phir, colat = checkinputs(lat, mlt)

    if (inside == 0):
        return fill

    phim: List[float] = [phir, phir * 2]
    cospm: List[float] = [cos(phi) for phi in phim]
    sinpm: List[float] = [sin(phi) for phi in phim]

    z = 0.
    skip: bool = False
    for j in range(reader.csize):
        if skip:
            skip = False
            continue

        m = reader.ms[j]
        if (reader.ab[j] == 1):
            plm = scplm(j, colat, nlm)
            skip = False
            if (m == 0):
                z = z + plm * esphc[j]
            else:
                z = z + plm * (esphc[j] * cospm[m] + esphc[j + 1] * sinpm[m])
                skip = True

    return z


def setboundary(angle: float, bt: float, tilt: float, swvel: float, swden: float, file_path: str):
    global tmat, ttmat, bndyfitr
    i: int

    swp: float
    # xc: float
    theta: float
    ct: float
    st: float
    tilt2: float
    cosa: float
    btx: float
    # x = np.zeros(reader.na)
    # c = np.zeros(reader.na)

    reader.read_bndy(file_path + '//W05scBndy.dat')

    # Calculate the transformation matrix to the coordinate system
    # of the offset pole.
    xc: float = 4.2
    theta = radians(xc)
    ct = cos(theta)
    st = sin(theta)

    tmat = [[ct, 0., st],
            [0., 1., 0.],
            [-st, 0., ct],
            ]
    ttmat = [[ct, 0., -st],
             [0., 1., 0.],
             [st, 0., ct],
             ]

    swp: float = swden * swvel ** 2 * 1.6726e-6  # pressure
    tilt2 = tilt ** 2
    cosa = cos(radians(angle))
    btx = 1. - exp(-bt * reader.ex_bndy[0])
    if (bt > 1.):
        btx = btx * bt ** reader.ex_bndy[1]
    else:
        cosa = 1. + bt * (cosa - 1.)  # remove angle dependency for IMF under 1 nT
    # endif
    x = [1., cosa, btx, btx * cosa, swvel, swp]
    c = reader.bndya
    bndyfitr = 0.
    for i in range(reader.na):
        bndyfitr = bndyfitr + x[i] * c[i]
    # enddo

    return


def setmodel(by: float, bz: float, tilt: float, swvel: float, swden: float, file_path: str, model: str):
    i: int
    j: int

    bt: float
    angle: float
    pi: float
    stilt: float
    stilt2: float
    sw: float
    swp: float
    swe: float
    c0: float
    rang: float
    cosa: float
    sina: float
    cos2a: float
    sin2a: float

    cfits: np.ndarray = np.zeros(reader.d1_pot, reader.csize)
    a = np.zeros(reader.d1_pot)

    assert model.strip() in ['epot', 'bpot'], "('>>> setmodel: model must be either epot or bpot')"

    if (model.strip() == 'epot'):
        reader.read_potential(file_path + '//W05scEpot.dat')
    else:
        reader.read_potential(file_path + '//W05scBpot.dat')

    reader.read_schatable(file_path + '//SCHAtable.dat')

    bt = sqrt(by ** 2 + bz ** 2)
    angle = degrees(atan2(by, bz))
    setboundary(angle, bt, tilt, swvel, swden, file_path)

    stilt = sin(degrees(tilt))
    stilt2 = stilt ** 2
    sw = bt * swvel / 1000.
    swe = (1. - exp(-sw * reader.ex_pot[1])) * sw ** reader.ex_pot[0]
    c0 = 1.
    swp = swvel ** 2 * swden * 1.6726e-6
    rang = radians(angle)
    cosa = cos(rang)
    sina = sin(rang)
    cos2a = cos(2. * rang)
    sin2a = sin(2. * rang)

    if (bt < 1.):  # remove angle dependency for IMF under 1 nT
        cosa = -1. + bt * (cosa + 1.)
        cos2a = 1. + bt * (cos2a - 1.)
        sina = bt * sina
        sin2a = bt * sin2a
        pass

    cfits = reader.schfits  # ! schfits(d1_pot, csize) is in module read_data

    a = [c0, swe, stilt, stilt2, swp,
         swe * cosa, stilt * cosa, stilt2 * cosa, swp * cosa,
         swe * sina, stilt * sina, stilt2 * sina, swp * sina,
         swe * cos2a, swe * sin2a]

    if (model.strip() == 'epot'):
        # esphc(:) = 0.
        for j in range(reader.csize):
            for i in range(reader.d1_pot):
                esphc[j] += cfits[i][j] * a[i]
            pass
        pass
        # write(6,"('setmodel: esphc=',/,(6e12.4))") esphc
    else:
        # bsphc(:) = 0.
        for j in range(reader.csize):
            for i in range(reader.d1_pot):
                bsphc[j] += cfits[i][j] * a[i]
                pass
            pass
        # pass
        # pass
        # write(6,"('setmodel: bsphc=',/,(6e12.4))") bsphc
    return
