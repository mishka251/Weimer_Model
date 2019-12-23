from typing import List
import numpy as np


def read_data():
    # Data read from W05scEpot.dat or W05scBpot.dat:
    csize: int = 28  # const
    d1_pot: int = 15  # const
    d2_pot: int = 18  # const

    ab: List[int] = np.zeros(csize)
    ls: List[int] = np.zeros(csize)
    ms: List[int] = np.zeros(csize)

    alschfits: List[List[float]] = np.zeros(d2_pot, csize)
    schfits: List[List[float]] = np.zeros(d1_pot, csize)
    ex_pot: List[List[float]] = np.zeros(2)

    maxl_pot: int
    maxm_pot: int

    # Data read from SCHAtable.dat:

    d1_scha: int = 19  # const
    d2_scha: int = 7  # const
    d3_scha: int = 68  # const

    allnkm: List[List[List[float]]] = np.zeros(d1_scha, d2_scha, d3_scha)

    maxk_scha: int
    maxm_scha: int

    th0s: List[float] = np.zeros(d3_scha)

    # Data read from W05scBndy.dat:

    na: int = 6  # const
    nb: int = 7  # const

    bndya: List[float] = np.zeros(na)
    bndyb: List[float] = np.zeros(nb)
    ex_bndy: List[float] = np.zeros(2)

    def read_potential(infile: str):
        """
        ! Read ascii data file W05scEpot.dat or W05scBpot.dat, written by
        !   pro write_potential (write_data.pro)
        """
        fname: str
        # Local:
        i: int
        lu: int = 20

        csize_rd: int
        d1_rd: int
        d2_rd: int

        file = open(infile)

        fname = file.read()
        ab: List[int] = [int(file.read()) for i in range(28)]

        csize_rd: int = int(file.read())
        d1_rd: int = int(file.read())
        d2_rd: int = int(file.read())

        assert csize_rd == csize, f"('read_potential: file '{infile}': incompatable csize: 'csize_rd=',{csize_rd},' csize=',{csize})"

        assert d1_rd == d1_pot, f"(' read_potential: file '{infile}': incompatable d1: 'd1_rd=',{d1_rd},' d1_pot='{d1_pot})"

        assert d2_rd == d2_pot, f"(' read_potential: file '{infile}': incompatable d2: ','d2_rd='{d2_rd}' d2_pot='{d2_pot})"

        for i in range(1, csize):
            for j in range(6):
                alschfits[i][j] = float(file.read())
                pass
            pass

        ex_pot = [float(file.read()), float(file.read())]

        for i in range(28):
            ls[i] = int(file.read())

        maxl_pot = int(file.read())
        maxm_pot = int(file.read())

        for i in range(28):
            ms[i] = int(file.read())

        for i in range(csize):
            for j in range(6):
                schfits[i][j] = float(file.read())
                pass
            pass

        file.close()
        return

    def read_schatable(infile: str):
        i: int
        j: int
        fname: str

        file = open(infile)
        fname = file.read()

        maxk_scha = int(file.read())
        maxm_scha = int(file.read())

        for i in range(d3_scha):
            for j in range(d2_scha):
                for k in range(6):
                    allnkm[i][j][k] = float(file.read())
                    pass
                pass
            pass

        for i in range(8):
            th0s[i] = float(file.read())

        return

    def read_bndy(infile:str):
        rd_na:int
        rd_nb:int

        file = open(infile)

        fname:str = file.read()
        rd_na = int(file.read())
        rd_nb = int(file.read())

        assert rd_na == na, f"('read_potential: file '{file}': incompatable na: ','rd_na='{rd_na}' na=',{na})"

        assert rd_nb == nb, f"('>>> read_potential: file '{file}': incompatable nb: ','rd_nb='{rd_nb}' nb={nb})"

        for i in range(8):
            bndya[i] = float(file.read())

        for i in range(8):
            bndyb[i] = float(file.read())

        for i in range(8):
            ex_bndy[i] = float(file.read())

        file.close()

