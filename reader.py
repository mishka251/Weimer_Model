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

    alschfits: List[List[float]] = []
    schfits: List[List[float]] = []
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

        fname: str = file.readline()
        line = file.readline()
        self.ab: List[int] = [int(digit) for digit in line.strip().split(' ') if digit != '']

        line = file.readline()
        [csize_rd, d1_rd, d2_rd] = [int(digit) for digit in line.strip().split(' ')]
        # csize_rd: int = int(file.read())
        # d1_rd: int = int(file.read())
        # d2_rd: int = int(file.read())

        assert csize_rd == self.csize, \
            f"('read_potential: file '{infile}': incompatable csize: 'csize_rd=',{csize_rd},' csize=',{self.csize})"

        assert d1_rd == self.d1_pot, \
            f"(' read_potential: file '{infile}': incompatable d1: 'd1_rd=',{d1_rd},' d1_pot='{self.d1_pot})"

        assert d2_rd == self.d2_pot, \
            f"(' read_potential: file '{infile}': incompatable d2: ','d2_rd='{d2_rd}' d2_pot='{self.d2_pot})"

        self.alschfits = []
        for i in range(self.d2_pot):
            self.alschfits.append([])

        col = []

        while len(self.alschfits[0]) < self.csize:  # for i in range(self.csize):

            elements = [float(digit) for digit in file.readline().strip().split(' ') if digit != '']
            col.extend(elements)

            if len(col) == self.d2_pot:
                for i in range(self.d2_pot):
                    self.alschfits[i].append(col[i])
                col = []

            # pass
            pass

        line = file.readline()
        self.ex_pot = [float(digit) for digit in line.strip().split(' ') if digit != '']

        self.ls = []
        while len(self.ls) < self.csize:
            self.ls.extend([int(char) for char in file.readline().strip().split() if char != ''])
        # for i in range(self.csize):
        # self.ls[i] = int(file.read())

        [self.maxl_pot, self.maxm_pot] = self.ls = [int(l) for l in file.readline().strip().split() if l != '']

        assert len(self.ls) == self.csize, "ls!=csize"
        assert len(self.ms) == self.csize, "ms!=csize"
        # for i in range(self.csize):
        #     self.ms[i] = int(file.read())
        self.schfits = []
        for i in range(self.d2_pot):
            self.schfits.append([])

        col = []

        while len(self.schfits[0]) < self.csize:  # for i in range(self.csize):

            elements = [float(digit) for digit in file.readline().strip().split(' ') if digit != '']
            col.extend(elements)

            if len(col) == self.d1_pot:
                for i in range(self.d1_pot):
                    self.schfits[i].append(col[i])
                col = []

            # pass
            pass

        file.close()
        print("f1 read")
        return

    def read_schatable(self, infile: str):
        file = open(infile)
        _: str = file.readline()

        [self.maxk_scha, self.maxm_scha] = [int(char) for char in file.readline().strip().split() if char != '']
        # = int(file.read())

        for i in range(self.d3_scha):
            for j in range(self.d2_scha):
                col = []
                while len(col) < self.d1_scha:
                    col.extend([float(char) for char in file.readline().strip().split() if char != ''])
                for k in range(self.d1_scha):
                    self.allnkm[k][j][i] = col[k]
                    pass
                pass
            pass

        assert len(self.allnkm) == self.d1_scha, "d1_scha"
        assert all(map(lambda l: len(l) == self.d2_scha, self.allnkm)), "d2_scha"
        assert all(map(lambda l2: all(map(lambda l: len(l) == self.d3_scha, l2)), self.allnkm)), "d3_scha"

        self.th0s = []
        while len(self.th0s) < self.d3_scha:
            self.th0s.extend([float(char) for char in file.readline().strip().split() if char != ''])
            # for i in range(8):
            #     self.th0s[i] = float(file.read())
        assert len(self.th0s) == self.d3_scha, "th0s"
        return

    def read_bndy(self, infile: str):
        file = open(infile)

        _: str = file.readline()
        [rd_na, rd_nb] = [int(char) for char in file.readline().strip().split() if char != '']

        assert rd_na == self.na, \
            f"('read_potential: file '{file}': incompatable na: ','rd_na='{rd_na}' na=',{self.na})"

        assert rd_nb == self.nb, \
            f"('>>> read_potential: file '{file}': incompatable nb: ','rd_nb='{rd_nb}' nb={self.nb})"

        self.bndya = [float(char) for char in file.readline().strip().split() if char != '']
        # for i in range(8):
        #     self.bndya[i] = float(file.read())
        assert len(self.bndya) == self.na

        self.bndyb = [float(char) for char in file.readline().strip().split() if char != '']

        assert len(self.bndyb) == self.nb
        # for i in range(8):
        #     self.bndyb[i] = float(file.read())
        self.ex_bndy = [float(char) for char in file.readline().strip().split() if char != '']
        assert len(self.ex_bndy) == 2
        # for i in range(8):
        #     self.ex_bndy[i] = float(file.read())

        file.close()
        return
