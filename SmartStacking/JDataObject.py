# Data Handling Module specifically DIAS to start
import numpy as np
# define variables here


class JinDipole:
    # object containing current information
    # option 1 create from raw
    # option 2 create from file

    def __init__(self, InDpInfo):
        try:
            self.Tx1 = float(InDpInfo.Tx1)
        except:
            pass
        self.Tx1East = float(InDpInfo.Tx1East)
        self.Tx1North = float(InDpInfo.Tx1North)
        self.Tx1Elev = float(InDpInfo.Tx1Elev)
        self.Tx2East = float(InDpInfo.Tx2East)
        self.Tx2North = float(InDpInfo.Tx2North)
        self.In = float(InDpInfo.In)
        self.In_err = float(InDpInfo.In_err)


class JvoltDipole:
    # object containing voltage information
    # option 1 create from raw
    # option 2 create from file
    def __init__(self, VoltDpinfo):
        self.dipole = VoltDpinfo.DIPOLE
        try:
            self.Rx1 = float(VoltDpinfo.Rx1)
        except:
            pass
        self.Rx1File = VoltDpinfo.Rx1File
        # self.Rx1Relay = VoltDpinfo.Rx1Relay
        self.Rx1East = float(VoltDpinfo.Rx1East)
        self.Rx1North = float(VoltDpinfo.Rx1North)
        self.Rx1Elev = float(VoltDpinfo.Rx1Elev)
        self.Rx2File = VoltDpinfo.Rx2File
        self.Rx2East = float(VoltDpinfo.Rx2East)
        self.Rx2North = float(VoltDpinfo.Rx2North)
        self.Rx2Elev = float(VoltDpinfo.Rx2Elev)
        # self.Rx2Relay = VoltDpinfo.Rx2Relay
        self.Vp = float(VoltDpinfo.Vp)
        self.Vp_err = float(VoltDpinfo.Vp_err)
        self.Rho = float(VoltDpinfo.Rho)
        self.flagRho = VoltDpinfo.Rho_QC
        self.Stack = float(VoltDpinfo.Stack)
        self.Mx = float(VoltDpinfo.Mx)
        self.Mx_err = float(VoltDpinfo.Mx_err)
        self.flagMx = VoltDpinfo.Mx_QC
        self.flagBad = VoltDpinfo.Status
        self.TimeBase = VoltDpinfo.TimeBase

        # somehow get all the decay info
        # TODO: must have smarts to tell how many windows
        # do this by using the length of the header ID's
        # Find where = M01 then subrtact from length
        self.Vs = np.zeros(20)
        self.Vs[0] = VoltDpinfo.Vs01
        self.Vs[1] = VoltDpinfo.Vs02
        self.Vs[2] = VoltDpinfo.Vs03
        self.Vs[3] = VoltDpinfo.Vs04
        self.Vs[4] = VoltDpinfo.Vs05
        self.Vs[5] = VoltDpinfo.Vs06
        self.Vs[6] = VoltDpinfo.Vs07
        self.Vs[7] = VoltDpinfo.Vs08
        self.Vs[8] = VoltDpinfo.Vs09
        self.Vs[9] = VoltDpinfo.Vs10
        self.Vs[10] = VoltDpinfo.Vs11
        self.Vs[11] = VoltDpinfo.Vs12
        self.Vs[12] = VoltDpinfo.Vs13
        self.Vs[13] = VoltDpinfo.Vs14
        self.Vs[14] = VoltDpinfo.Vs15
        self.Vs[15] = VoltDpinfo.Vs16
        self.Vs[16] = VoltDpinfo.Vs17
        self.Vs[17] = VoltDpinfo.Vs18
        self.Vs[18] = VoltDpinfo.Vs19
        self.Vs[19] = VoltDpinfo.Vs20

    def getXplotpoint(self, Idp):
        if (self.Rx1 > Idp.Tx1):
            x = Idp.Tx1 + ((self.Rx1 - Idp.Tx1) / 2.0)
        elif (self.Rx1 < self.Idp):
            x = Idp.Tx1 - ((Idp.Tx1 - self.Rx1) / 2.0)
        return[x]

    def getZplotpoint(self, Idp):
        z = -(abs(Idp.Tx1 - self.Rx1)) / 2.0
        return[z]

    def calcRho(self, Idp):
        arho = 10
        return[arho]


class Jreading:
    # object to handle current and voltage dipole
    # information for a given reading
    def __init__(self, mem):
        self.MemNumber = mem
        self.Vdp = []
    # method for adding voltage dipoles

    def addVoltageDipole(self, JVdp):
        self.Vdp.append(JVdp)

    # method for assigning Current dipole
    def addInDipole(self, JIdp):
        self.Idp = JIdp


class Jpatch:
    # object to handle reading for a given data patch
    def __init__(self, Id):
        self.ID = Id
        self.readings = []

    def addreading(self, Jrdg):
        self.readings.append(Jrdg)


class Jreadtxtline:
    def __init__(self, hdrLine, dataLine):
        # make a structure with the header inputs
        for n in range(len(hdrLine)):
            setattr(self, hdrLine[n], dataLine[n])
