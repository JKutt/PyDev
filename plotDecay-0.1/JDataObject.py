#################################################
# Imports

import numpy as np
# from scipy.special import factorial
# from SimPEG.EM.Static import DC

##################################################
# define methods


def loadDias(fileName):
    """
    Function for loading a dias file and returns a
    "Patch" class complete with sources and recievers

    Input:
    fileName = complete path to data file

    """
    lines = 0
    text_file = open(fileName, "r")

    # determin how many lines in the file
    while text_file.readline():
            lines += 1
    text_file.close()

    # initiate a patch
    patch = Jpatch()
    # read header information
    text_file = open(fileName, "r")
    # initiate reading control variables
    currRdg = 0
    for i, line in enumerate(text_file):
        if i == 4:
            Varinfo = line.split()
            header4 = line
            # print(Varinfo)
        elif i == 0:
                header1 = line
        elif i == 1:
                header2 = line
                id_info = line.split()
                patch.assignPatchID(id_info[1])
        elif i == 2:
                header3 = line
        elif i > 3:
            try:
                    datatxt = line.split()
                    # do some Jdatamanagment stuff
                    # print(datatxt)
                    varFields = Jreadtxtline(Varinfo, datatxt)
                    # verify if line is a new reading
                    if varFields.RDG == currRdg:
                        # add the dipoles
                        # Idp = Jdata.JinDipole(varFields)
                        Vpdp = JvoltDipole(varFields)
                        Rdg.addVoltageDipole(Vpdp)
                    else:
                        # create a reading
                        Rdg = Jreading(varFields.RDG)
                        Idp = JinDipole(varFields)
                        Vpdp = JvoltDipole(varFields)
                        Rdg.addVoltageDipole(Vpdp)
                        Rdg.addInDipole(Idp)
                        # add reading to the patch
                        patch.addreading(Rdg)
                        currRdg = varFields.RDG
            except:
                    pass

    text_file.close()
    headers = [header1, header2, header3, header4]
    patch.assignHeaderInfo(headers)
    return patch

######################################################
# Define Classes


class dcipTimeSeries:
    """
    Class containing Source information

    Input:
    sample_rate = number of samples/sec
    time_base = period of base frequency
    data = a numpy array containing the data
    """
    def __init__(self, sample_rate, time_base, data):
        self.timebase = time_base
        self.samplerate = sample_rate
        self.data
        # determine samples per period and half period
        self.sampPerT = time_base * sample_rate
        self.sampPerHalfT = self.sampPerT / 2.0

    def stack(self):
        """
        simple brute stack algorithm

        """

    def stackTimeSeries(self, filterKernal):
        """
        TODO
        Input:
        filterKernal = a numpy array consisting of filter kernal
        e.g Hanning, Kaiser, etc...
        """
        if filterKernal.size > 3:
            bkernal = np.ones((3, 1))              # create basic 3 T kernal
            bkernal[1] = -2.0
            bsgn = np.ones((1, filterKernal.size))
            bsgn[0, 0::2] = bsgn[0, 0::2] * -1.0         # linear drift removal
            bwt = np.zeros((3, filterKernal.size))
            bwt = np.matmul(bkernal, bsgn * filterKernal) # creates stack kernal

            # map and sum the weighted kernal of bwt
            k0 = np.range(filterKernal.size * 3) + 1
            tmp1 = np.arange(3) + 1
            tmp2 = np.ones(filterKernal.size)
            tmp3 = np.ones((3, 1))
            tmp4 = np.arange(filterKernal.size)
            knew = np.matmul(tmp1, tmp2) + np.matmul(tmp3, (tmp4 * (filterKernal.size + 3)))

    def ensembleStackTimeSeries(self, filterKernal):
        """
        TODO
        Input:
        filterKernal = a numpy array consisting of filter kernal
                note: the size of the filter Kernal will be used
                to determine how many ensembles.
                e.g Hanning, Kaiser, etc...
        error_allowance = std of acceptance of decay
        """
        self.data = 100. / 1000

    def statRejectTimeSeries(self, error_allowance):
        """
        TODO
        Input:
        filterKernal = a numpy array consisting of filter kernal
                e.g Hanning, Kaiser, etc...
        error_allowance = std of acceptance of decay
        """
        self.data = 100. / 1000


# ===================================================
# Dias Data specific class
class JinDipole:
    """
    Class containing Source information

    Initiate with a structure containing location
    and Current value + error

    """

    def __init__(self, InDpInfo):
        try:
            self.Tx1 = float(InDpInfo.Tx1x)
        except:
            pass
        self.Tx1East = float(InDpInfo.Tx1East)
        self.Tx1North = float(InDpInfo.Tx1North)
        self.Tx1Elev = float(InDpInfo.Tx1Elev)
        self.Tx2East = float(InDpInfo.Tx2East)
        self.Tx2North = float(InDpInfo.Tx2North)
        self.Tx2Elev = float(InDpInfo.Tx2Elev)
        self.In = float(InDpInfo.In)
        self.In_err = float(InDpInfo.In_err)


class JvoltDipole:
    """
    object containing voltage information

    """

    def __init__(self, VoltDpinfo):
        self.dipole = VoltDpinfo.DIPOLE
        try:
            self.Rx1 = float(VoltDpinfo.Rx1x)
        except:
            pass
        try:
            self.Rx2 = float(VoltDpinfo.Rx2x)
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
        self.couple = float(VoltDpinfo.Coupling)
        self.Stack = float(VoltDpinfo.Stack)
        try:
            self.Mx = float(VoltDpinfo.Mx)
        except:
            self.Mx = -99.9
        self.Mx_err = float(VoltDpinfo.Mx_err)
        self.flagMx = VoltDpinfo.Mx_QC
        self.flagBad = VoltDpinfo.Status
        self.TimeBase = VoltDpinfo.TimeBase
        self.Vs = np.asarray(VoltDpinfo.Vs)

    def getXplotpoint(self, Idp):
        if (self.Rx1 > Idp.Tx1):
            x = Idp.Tx1 + ((self.Rx1 - Idp.Tx1) / 2.0)
        elif (self.Rx1 < Idp.Tx1):
            x = Idp.Tx1 - ((Idp.Tx1 - self.Rx1) / 2.0)
        else:
            print
        return[x]

    def getZplotpoint(self, Idp):
        z = -(abs(Idp.Tx1 - self.Rx1)) / 2.0
        return[z]

    def calcRho(self):
        return[self.Rho]


class Jreading:
    """
    Class to handle current and voltage dipole
    information for a given source

    """
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
    """
    Class to hold source information for a given data patch

    """

    def __init__(self):
        self.readings = []

    def addreading(self, Jrdg):
        self.readings.append(Jrdg)

    def assignPatchID(self, id):
        self.ID = id

    def assignHeaderInfo(self, headerLines):
        """
        Method for processing the header lines
        of a Dias data file. (e.g IP times, project name, etc.)

        Input: an array of header lines from file

        """

        self.headers = headerLines         # assigns headers to patch class
        # process IP times from header
        timing_string = self.headers[2].split(' ')[1]
        timings = timing_string.split(';')[0]
        timing = timings.split(',')
        self.window_start = np.zeros(len(timing))
        self.window_end = np.zeros(len(timing))
        self.window_center = np.zeros(len(timing))
        self.window_width = np.zeros(len(timing))
        # loops each available IP window and calculate width and centre
        for i in range(len(timing)):
            wintime = timing[i].split(':')
            self.window_start[i] = float(wintime[0])
            self.window_end[i] = float(wintime[1])
            self.window_center[i] = (self.window_start[i] +
                                     self.window_end[i]) / 2.0
            self.window_width[i] = (self.window_end[i] -
                                    self.window_start[i])

    # def createDcSurvey(self, data_type):
    #     """
    #     Loads a dias data file to a SimPEG "srcList" class

    #     Input:
    #     datatype = Choose either IP or DC

    #     Note: elevation is -ve for simPEG inversion

    #     """

    #     srcLists = []                                 # Injections + dipoles
    #     data = []                                     # data from file
    #     d_weights = []                                # weights for the data
    #     xpp = []
    #     ypp = []
    #     num_rdg = len(self.readings)
    #     minE = self.readings[0].Idp.Tx2East
    #     minN = self.readings[0].Idp.Tx2North
    #     maxN = minN
    #     maxE = minE
    #     for k in range(num_rdg):
    #         num_dipole = len(self.readings[k].Vdp)
    #         rx = np.zeros((num_dipole, 6))
    #         tx = np.array([self.readings[k].Idp.Tx1East,
    #                       self.readings[k].Idp.Tx1North,
    #                       -self.readings[k].Idp.Tx1Elev,
    #                       self.readings[k].Idp.Tx2East,
    #                       self.readings[k].Idp.Tx2North,
    #                       -self.readings[k].Idp.Tx2Elev])
    #         if self.readings[k].Idp.Tx1East > maxE:
    #             maxE = self.readings[k].Idp.Tx1East
    #         if self.readings[k].Idp.Tx1East < minE:
    #             minE = self.readings[k].Idp.Tx1East
    #         if self.readings[k].Idp.Tx1North > maxN:
    #             maxN = self.readings[k].Idp.Tx1North
    #         if self.readings[k].Idp.Tx1North < minN:
    #             minN = self.readings[k].Idp.Tx1North

    #         for i in range(num_dipole):
    #             try:
    #                 xpp.append(self.readings[k].Idp.Tx1East)
    #                 ypp.append(self.readings[k].Idp.Tx1Elev)
    #                 rx[i, :] = [self.readings[k].Vdp[i].Rx1East,
    #                             self.readings[k].Vdp[i].Rx1North,
    #                             -self.readings[k].Vdp[i].Rx1Elev,
    #                             self.readings[k].Vdp[i].Rx2East,
    #                             self.readings[k].Vdp[i].Rx2North,
    #                             -self.readings[k].Vdp[i].Rx2Elev]
    #                 if data_type == "DC":
    #                     data.append(self.readings[k].Vdp[i].Vp /
    #                                 self.readings[k].Idp.In)
    #                     d_weights.append((self.readings[k].Vdp[i].Vp /
    #                                      self.readings[k].Idp.In) *
    #                                      (self.readings[k].Vdp[i].Vp_err /
    #                                       100.0))
    #                 if data_type == "IP":
    #                     data.append(self.readings[k].Vdp[i].Vp /
    #                                 self.readings[k].Idp.In)
    #                     d_weights.append((self.readings[k].Vdp[i].Mx *
    #                                      (self.readings[k].Vdp[i].Mx_err /
    #                                       100.0)))
    #             except:
    #                 pass

    #         Rx = DC.Rx.Dipole(rx[:, :3], rx[:, 3:])    # create dipole list
    #         srcLists.append(DC.Src.Dipole([Rx], tx[:3], tx[3:]))

    #     survey = DC.SurveyDC.Survey(srcLists)          # creates the survey
    #     survey.dobs = np.asarray(data)                 # assigns data
    #     survey.std = np.asarray(d_weights)             # assign data weights
    #     survey.eps = 0.

    #     return {'dc_survey': survey}


class Jreadtxtline:
    """
    Class specifically for reading a line of text from a dias file

    """

    def __init__(self, hdrLine, dataLine):
        # make a structure with the header inputs
        self.Vs = []
        for n in range(len(hdrLine)):
            if hdrLine[n].find("Vs") == 0:
                self.Vs.append(float(dataLine[n]))
            else:
                setattr(self, hdrLine[n], dataLine[n])
