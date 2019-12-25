from math import sqrt, atan2, radians, degrees, sin, exp, cos, atan, log, asin, factorial
from typing import List, Tuple

from reader import Reader
import numpy as np
from utils import value_locate, interpol_quad, lngamma, km_n


class Calculator:
    reader: Reader = Reader()

    bndyfitr: float
    esphc: List[float] = [0] * reader.csize
    bsphc: List[float] = [0] * reader.csize
    tmat: List[List[float]] = np.zeros((3, 3), np.float)
    ttmat: List[List[float]] = np.zeros((3, 3), np.float)
    #
    mxtablesize: int = 200

    #
    plmtable: List[List[float]] = np.zeros((mxtablesize, reader.csize))
    colattable: List[float] = np.zeros(mxtablesize, np.float)
    #
    nlms: List[float] = np.zeros(reader.csize, np.float)

    def do_rotation(self, latin: float, lonin: float) -> Tuple[float, float]:
        """
        :param latin:
        :param lonin:
        :return:latout,lonout
        """
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
            pos[i] = self.tmat[0][i] * a + self.tmat[1][i] * b + self.tmat[2][i] * stc

        latout: float = degrees(asin(pos[0]))
        lonout: float = degrees(atan2(pos[1], pos[0]))

        return latout, lonout

    def checkinputs(self, lat: float, mlt: float) -> Tuple[int, float, float]:
        """

        :param lat:
        :param mlt:
        :return:  inside: int, phir: float, colat: float
        """
        tlat: float
        tlon: float

        lon: float = mlt * 15.
        tlat, tlon = self.do_rotation(lat, lon)
        radii: float = 90. - tlat

        inside: int = 1 if (radii <= self.bndyfitr) else 0
        phir: float = radians(tlon)
        colat: float = radii

        return inside, phir, colat

    def nkmlookup(self, k: int, m: int, th0: float) -> float:
        if th0 == 90.:
            return float(k)

        kk: int = k + 1
        mm: int = m + 1

        th0a: float = th0

        res: List[float] = []
        max_k = self.reader.maxk_scha
        max_m = self.reader.maxm_scha
        th0s = self.reader.th0s
        all_nkm = self.reader.allnkm

        if kk >= max_k:
            print(f"('>>> nkmlookup: kk > maxk: kk='{kk}' maxk='{max_k}")
            res = interpol_quad(all_nkm[max_k - 1][mm], th0s, [th0a])
        if mm >= max_m:
            print(f"('>>> nkmlookup: mm > maxm: kk='{kk}' maxm='{max_m}")
            res = interpol_quad(all_nkm[kk][max_m - 1], th0s, [th0a])
        if th0 < th0s[0]:
            print(f"('>>> nkmlookup: th0 < th0s(1): th0='{th0}' th0s(1)='{th0s[0]}")

        res = interpol_quad(all_nkm[kk][mm], th0s, [th0a])
        return res[0]

    def pm_n(self, m: int, r: float, cth: List[float], table_size: int) -> List[float]:
        a: List[float]
        if m == 0:
            a = [0]
        else:
            a = [sqrt(1 - cth[i] ** 2) ** m for i in range(1, table_size + 1)]
        xn: float = r * (r + 1)

        x = [(1 - ct) / 2 for ct in cth]

        table: List[float] = a

        tmp: List[float] = [10000] * table_size
        k = 0
        while max(tmp) > 1e-6:
            for i in range(table_size):
                a[i] *= (x[i] * ((k + m - 1.) * (k + m) - xn) / (k * (k + m)))
                table[i] += a[i]
                pass

            k += 1

            for i in range(table_size):
                div = abs(table[i])
                div = max(div, 1e-6)
                tmp[i] = abs(a[i]) / div
                pass
            pass

        ans = km_n(m, r)
        return [t * ans for t in table]

    prev_th0 = 1e36

    def scplm(self, index: int, colat: float) -> Tuple[float, float]:
        ls = self.reader.ls
        ms = self.reader.ms
        ab = self.reader.ab
        # nlms:List[float] = [0]
        skip: bool = False
        th0 = self.bndyfitr
        if self.prev_th0 != th0:
            self.tablesize = 3 * int(round((th0)))
            assert self.tablesize <= self.mxtablesize, \
                f"('>>> tablesize > mxtablesize: tablesize={self.tablesize} mxtablesize={self.mxtablesize} tn0={th0}"

            colattable = [float(i - 1) * (th0 / float(self.tablesize - 1)) for i in range(self.tablesize)]
            cth = [cos(degrees(col)) for col in colattable]
            self.prev_th0 = th0
            nlms = [0] * self.reader.csize

            for j in range(self.reader.csize):
                if skip:
                    skip = False
                    continue
                    pass
                nlms.append(self.nkmlookup(ls[j], ms[j], th0))
                tmp: List[float] = self.pm_n(ms[j], nlms[j], cth, self.tablesize)
                for _i in range(self.tablesize):
                    self.plmtable[_i][index] = tmp[_i]
                    pass
                skip = False

                if ms[j] == 0 and ab[j] > 0:
                    self.plmtable[0][j + 1] = self.plmtable[0][j]
                    nlms[j + 1] = nlms[j]
                    skip = True
                    pass
                pass  # end for
            pass  # endif

        nlm = self.nlms[index]
        colata = [colat]

        tmp = [self.plmtable[i][index] for i in range(self.tablesize)]

        out = interpol_quad(tmp, self.colattable[1: self.tablesize], colata)
        return out[1], nlm

    def mpfac(self, lat: float, mlt: float, fill: float) -> float:
        """

        :param lat:
        :param mlt:
        :param fill:
        :return: fac:float
        """
        ls = self.reader.ls
        ms = self.reader.ms
        ab = self.reader.ab
        csize = self.reader.csize

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

        inside, phir, colat = self.checkinputs(lat, mlt)

        if (inside == 0):
            return fill

        phim: List[float] = [phir, phir * 2]
        cospm: List[float] = [cos(phi) for phi in phim]
        sinpm: List[float] = [sin(phi) for phi in phim]

        z = 0.
        # jloop: do
        skip: bool = False
        for j in range(csize):
            if skip:
                skip = False
                continue

            if ls[j] >= 11:
                break
            # m = ms[j]
            if ab[j] == 1:
                plm, nlm = self.scplm(j, colat)
                plm = plm * (nlm * (nlm + 1.))

                if ms[j] == 0:
                    z = z - plm * self.bsphc[j]
                else:
                    z = z - (plm * (self.bsphc[j] * cospm[ms[j]] + self.bsphc[j + 1] * sinpm[ms[j]]))
                    skip = True

        pi = 4. * atan(1.)
        cfactor = -1.e5 / (4. * pi * re ** 2)  # convert to uA/m2
        z = z * cfactor
        return z

    def epotval(self, lat: float, mlt: float, fill: float) -> float:
        """

        :param lat:
        :param mlt:
        :param fill:
        :param epot:
        :return: erot:float
        """

        # phir: float
        plm: float
        # colat: float
        nlm: float

        csize = self.reader.csize
        ms = self.reader.ms
        ab = self.reader.ab
        inside, phir, colat = self.checkinputs(lat, mlt)

        if (inside == 0):
            return fill

        phim: List[float] = [phir, phir * 2]
        cospm: List[float] = [cos(phi) for phi in phim]
        sinpm: List[float] = [sin(phi) for phi in phim]

        z: float = 0.
        skip: bool = False
        for j in range(csize):
            if skip:
                skip = False
                continue

            m: int = ms[j]
            if ab[j] == 1:
                plm, nlm = self.scplm(j, colat)
                skip = False
                if m == 0:
                    z = z + plm * self.esphc[j]
                else:
                    z = z + plm * (self.esphc[j] * cospm[m] + self.esphc[j + 1] * sinpm[m])
                    skip = True

        return z

    def setboundary(self, angle: float, bt: float, tilt: float, swvel: float, swden: float, file_path: str):
        btx: float
        self.reader.read_bndy(file_path + '//W05scBndy.dat')

        # Calculate the transformation matrix to the coordinate system
        # of the offset pole.
        xc: float = 4.2
        theta: float = radians(xc)
        ct: float = cos(theta)
        st: float = sin(theta)

        self.tmat = [[ct, 0., st],
                     [0., 1., 0.],
                     [-st, 0., ct],
                     ]

        self.ttmat = [[ct, 0., -st],
                      [0., 1., 0.],
                      [st, 0., ct],
                      ]

        swp: float = swden * swvel ** 2 * 1.6726e-6  # pressure
        self.tilt2: float = tilt ** 2
        cosa: float = cos(radians(angle))
        btx: float = 1. - exp(-bt * self.reader.ex_bndy[0])
        if (bt > 1.):
            btx = btx * bt ** self.reader.ex_bndy[1]
        else:
            cosa = 1. + bt * (cosa - 1.)  # remove angle dependency for IMF under 1 nT
        # endif
        x = [1., cosa, btx, btx * cosa, swvel, swp]
        c = self.reader.bndya
        self.bndyfitr = sum([x_i * c_i for x_i, c_i in zip(x, c)])

        return

    def setmodel(self, by: float, bz: float, tilt: float, swvel: float, swden: float, file_path: str, model: str):

        assert model.strip() in ['epot', 'bpot'], "('>>> setmodel: model must be either epot or bpot')"

        if (model.strip() == 'epot'):
            self.reader.read_potential(file_path + '//W05scEpot.dat')
        else:
            self.reader.read_potential(file_path + '//W05scBpot.dat')

        self.reader.read_schatable(file_path + '//SCHAtable.dat')

        bt: float = sqrt(by ** 2 + bz ** 2)
        angle: float = degrees(atan2(by, bz))
        self.setboundary(angle, bt, tilt, swvel, swden, file_path)

        stilt: float = sin(degrees(tilt))
        stilt2: float = stilt ** 2
        sw: float = bt * swvel / 1000.
        swe: float = (1. - exp(-sw * self.reader.ex_pot[1])) * sw ** self.reader.ex_pot[0]
        c0: float = 1.
        swp: float = swvel ** 2 * swden * 1.6726e-6
        rang: float = radians(angle)
        cosa: float = cos(rang)
        sina: float = sin(rang)
        cos2a: float = cos(2. * rang)
        sin2a: float = sin(2. * rang)

        if bt < 1.:  # remove angle dependency for IMF under 1 nT
            cosa = -1. + bt * (cosa + 1.)
            cos2a = 1. + bt * (cos2a - 1.)
            sina = bt * sina
            sin2a = bt * sin2a
            pass

        cfits: List[List[float]] = self.reader.schfits  # ! schfits(d1_pot, csize) is in module read_data

        a: List[float] = [c0, swe, stilt, stilt2, swp,
                          swe * cosa, stilt * cosa, stilt2 * cosa, swp * cosa,
                          swe * sina, stilt * sina, stilt2 * sina, swp * sina,
                          swe * cos2a, swe * sin2a]

        if (model.strip() == 'epot'):
            # esphc(:) = 0.
            for j in range(self.reader.csize):
                for i in range(self.reader.d1_pot):
                    self.esphc[j] += cfits[i][j] * a[i]
                pass
            pass
            # write(6,"('setmodel: esphc=',/,(6e12.4))") esphc
        else:
            # bsphc(:) = 0.
            for j in range(self.reader.csize):
                for i in range(self.reader.d1_pot):
                    self.bsphc[j] += cfits[i][j] * a[i]
                    pass
                pass
            # pass
            # pass
            # write(6,"('setmodel: bsphc=',/,(6e12.4))") bsphc
        return
