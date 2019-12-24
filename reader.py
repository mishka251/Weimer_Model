from typing import List, TextIO
import numpy as np


# Data read from W05scEpot.dat or W05scBpot.dat:

class Reader:
    csize: int = 28  # const
    d1_pot: int = 15  # const
    d2_pot: int = 18  # const

    ab: List[int] = np.zeros(csize)
    ls: List[int] = np.zeros(csize)
    ms: List[int] = np.zeros(csize)

    alschfits: List[List[float]] = np.zeros((d2_pot, csize), np.float)
    schfits: List[List[float]] = np.zeros((d1_pot, csize), np.float)
    ex_pot: List[float] = np.zeros(2, np.float)

    maxl_pot: int
    maxm_pot: int

    # Data read from SCHAtable.dat:

    d1_scha: int = 19  # const
    d2_scha: int = 7  # const
    d3_scha: int = 68  # const

    allnkm: List[List[List[float]]] = np.zeros((d1_scha, d2_scha, d3_scha), np.float)

    maxk_scha: int
    maxm_scha: int

    th0s: List[float] = np.zeros(d3_scha)

    # Data read from W05scBndy.dat:

    na: int = 6  # const
    nb: int = 7  # const

    bndya: List[float] = np.zeros(na)
    bndyb: List[float] = np.zeros(nb)
    ex_bndy: List[float] = np.zeros(2)

    def read_potential(self, infile: str):
        """
         Read ascii data file W05scEpot.dat or W05scBpot.dat, written by
        pro write_potential (write_data.pro)
        """

        file: TextIO = open(infile)

        _: str = file.read()
        self.ab: List[int] = [int(file.read()) for i in range(28)]

        csize_rd: int = int(file.read())
        d1_rd: int = int(file.read())
        d2_rd: int = int(file.read())

        assert csize_rd == self.csize, \
            f"('read_potential: file '{infile}': incompatable csize: 'csize_rd=',{csize_rd},' csize=',{self.csize})"

        assert d1_rd == self.d1_pot, \
            f"(' read_potential: file '{infile}': incompatable d1: 'd1_rd=',{d1_rd},' d1_pot='{self.d1_pot})"

        assert d2_rd == self.d2_pot, \
            f"(' read_potential: file '{infile}': incompatable d2: ','d2_rd='{d2_rd}' d2_pot='{self.d2_pot})"

        for i in range(self.csize):
            for j in range(6):
                self.alschfits[i][j] = float(file.read())
                pass
            pass

        self.ex_pot = [float(file.read()), float(file.read())]

        for i in range(self.csize):
            self.ls[i] = int(file.read())

        self.maxl_pot = int(file.read())
        self.maxm_pot = int(file.read())

        for i in range(self.csize):
            self.ms[i] = int(file.read())

        for i in range(self.csize):
            for j in range(6):
                self.schfits[i][j] = float(file.read())
                pass
            pass

        file.close()
        return

    def read_schatable(self, infile: str):
        file = open(infile)
        _: str = file.read()

        self.maxk_scha = int(file.read())
        self.maxm_scha = int(file.read())

        for i in range(self.d3_scha):
            for j in range(self.d2_scha):
                for k in range(6):
                    self.allnkm[i][j][k] = float(file.read())
                    pass
                pass
            pass

        for i in range(8):
            self.th0s[i] = float(file.read())

        return

    def read_bndy(self, infile: str):
        file = open(infile)

        _: str = file.read()
        rd_na = int(file.read())
        rd_nb = int(file.read())

        assert rd_na == self.na, \
            f"('read_potential: file '{file}': incompatable na: ','rd_na='{rd_na}' na=',{self.na})"

        assert rd_nb == self.nb, \
            f"('>>> read_potential: file '{file}': incompatable nb: ','rd_nb='{rd_nb}' nb={self.nb})"

        for i in range(8):
            self.bndya[i] = float(file.read())

        for i in range(8):
            self.bndyb[i] = float(file.read())

        for i in range(8):
            self.ex_bndy[i] = float(file.read())

        file.close()
        return
