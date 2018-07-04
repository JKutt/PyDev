#################################################
# Imports
import numpy as np
##################################################
# define methods


def calcColeCole(mx_decay, window_widths, cond_init, tau_init, radius):
    """
    takes in a decay and returns best fitting cole-cole
    note:
    A tool to calculate the time - domain voltage response[V(t) / Vo] for a
    homogenous Cole - Cole model of the earth using the digital linear filter
    formulation given by Guptasarma(Geophys, vol 47, pg 1575, 1982)
    """
    # check polarity
    polarity_check = np.sum(mx_decay)
    if polarity_check < 0.0:
        mx_decay = mx_decay * -1.0

    time = np.zeros(window_widths.size)  # initiates time
    # convert window widths to accumlative times specific to algorithm
    for i in range(time.size):
        time[i] = (np.sum(window_widths[0:i + 1]) / 2.0) / 1000.0
    c = np.zeros(9)                      # conductivity array
    v_cole = np.zeros(window_widths.size)  # best fit cole-cole
    tau = np.zeros(9)                    # time constant array
    err = np.zeros((9, 9))               # error matrix
    cole_m = np.zeros((9, 9))            # matrix of chargeabilities
    radius = radius                        # radius of array fill
    c[5] = cond_init                          # center of cond. array
    tau[5] = tau_init                   # center of tau array
    tau10 = np.log10(tau[5])             # log of time constant
    idx = np.arange(0, 9)
    c = c[5] + radius * (idx - 5) / 40.0  # fill cond. array
    tau = np.power(10.0,
                   (tau10 + radius * (idx - 5) / 2.))  # fill tau array
    # create filter
    areg = np.asarray([-3.82704, -3.56608, -3.30512, -3.04416,
                       -2.78320, -2.52224, -2.26128, -2.00032,
                       -1.73936, -1.47840, -1.21744, -.95648,
                       -.69552, -0.43456, -0.17360, 0.08736,
                       0.34832, 0.60928, 0.87024, 1.13120, 1.39216])
    # create 2nd filter
    preg = np.asarray([0.000349998, -0.000418371, 0.000772828,
                      -0.000171356, 0.001022172, 0.000897638,
                      0.002208974, 0.003844944, 0.006809040,
                      0.013029162, 0.022661391, 0.042972904,
                      0.075423603, 0.139346367, 0.234486236,
                      0.366178323, 0.284615486, -0.235691746,
                      0.046994188, -0.005901946, 0.000570165])
    fit_weights = np.ones(time.size) + 0.5  # create filter weights
    fit_weights[0] = 0.3
    fit_weights[1] = 0.5
    fit_weights[2] = 0.5
    v_cole = np.zeros(time.size)      # initiate decay array
    minErr = 0.01                      # signify initial Low error
    c_idx = 0                         # index of cond. of min err
    tau_idx = 0                       # index of cond. of min err

    # loop through the arrays of cond. and tau
    for i in range(c.size):
        for j in range(tau.size):
            ax = c[5] * np.pi / 2.0
            for win in range(mx_decay.size):
                v_temp = 0.0
                for n in range(areg.size):
                    w = np.power(10.0, (areg[n] - np.log10(time[win])))
                    ex = np.power(w * tau[j], c[i])
                    y = np.complex(ex * np.cos(ax), ex * np.sin(ax))
                    z = 1.0 - 1.0 / (1.0 + y)
                    v_temp = v_temp + preg[n] * np.real(z)
                v_cole[win] = v_temp

            # calculate error
            # norm_weights = np.sum(fit_weights) / fit_weights.size
            # print(norm_weights)
            serr = (np.sum(np.power((mx_decay - v_cole), 2) *
                    fit_weights))
            err = np.sqrt(serr / (fit_weights.size))
            if err < minErr:
                c_idx = i
                tau_idx = j
                minErr = (err)
                # print(serr)
                # print(err)

            mx_conversion = 1000.                 # conversion for mV/V
            cole_m[i, j] = (np.sum(v_cole * window_widths) /
                            np.sum(window_widths)) * mx_conversion  # calcs Mx

            # go back and calculate best fit cole-cole curve and save it
            for win in range(mx_decay.size):
                v_temp = 0.0
                for n in range(areg.size):
                    w = np.power(10.0, (areg[n] - np.log10(time[win])))
                    ex = np.power(w * tau[tau_idx], c[c_idx])
                    y = np.complex(ex * np.cos(ax), ex * np.sin(ax))
                    z = 1.0 - 1.0 / (1.0 + y)
                    v_temp = v_temp + preg[n] * np.real(z)
                v_cole[win] = v_temp

            # calculate the percent diff
            # percent_diff = (
            #     np.mean((np.abs((v_cole - mx_decay) /
            #             ((v_cole + mx_decay) / 2.))) * fit_weights))

    return c[c_idx], tau[tau_idx], cole_m[c_idx, tau_idx], minErr / 3. * 1000.


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
        if i == 3:
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
    """
    object containing voltage information

    """

    def __init__(self, VoltDpinfo):
        self.RDG = VoltDpinfo.RDG
        self.dipole = VoltDpinfo.DIPOLE
        try:
            self.Rx1 = float(VoltDpinfo.Rx1)
        except:
            pass
        try:
            self.Rx2 = float(VoltDpinfo.Rx2)
        except:
            pass
        self.Rx1File = VoltDpinfo.Rx1FileName
        self.Status = VoltDpinfo.Status
        self.Nom_a = float(VoltDpinfo.Nominal_a)
        self.n = VoltDpinfo.n
        # self.Rx1Relay = VoltDpinfo.Rx1Relay
        self.Rx1East = float(VoltDpinfo.Rx1East)
        self.Rx1North = float(VoltDpinfo.Rx1North)
        self.Rx1Elev = float(VoltDpinfo.Rx1Elev)
        self.Rx2File = VoltDpinfo.Rx2FileName
        self.Rx2East = float(VoltDpinfo.Rx2East)
        self.Rx2North = float(VoltDpinfo.Rx2North)
        self.Rx2Elev = float(VoltDpinfo.Rx2Elev)
        # self.Rx2Relay = VoltDpinfo.Rx2Relay
        self.contact = float(0.00)
        self.Vp = float(VoltDpinfo.Vp)
        self.Vp_err = float(VoltDpinfo.Vp_err)
        self.Rho = float(VoltDpinfo.Rho)
        self.flagRho = VoltDpinfo.Rho_QC
        self.Stack = float(VoltDpinfo.Stack)
        try:
            self.Mx = float(VoltDpinfo.MA)
        except:
            self.Mx = -99.9
        self.Mx_err = float(VoltDpinfo.MA_err)
        self.flagMx = VoltDpinfo.MA_QC
        self.flagBad = VoltDpinfo.Status
        self.TimeBase = VoltDpinfo.TimeBase
        self.Vs = np.asarray(VoltDpinfo.Vs)
        self.cole_c = -99.9
        self.cole_tau = -99.9
        self.cole_m = -99.9

    def getXplotpoint(self, Idp):
        if (self.Rx1 > Idp.Tx1):
            x = Idp.Tx1 + ((self.Rx1 - Idp.Tx1) / 2.0)
        elif (self.Rx1 < Idp.Tx1):
            x = Idp.Tx1 - ((Idp.Tx1 - self.Rx1) / 2.0)
        return[x]

    def getZplotpoint(self, Idp):
        z = -(abs(Idp.Tx1 - self.Rx1)) / 2.0
        return[z]

    def calcRhoError(self, Idp):
        rho_error = 0.
        r1 = (np.power((np.power(Idp.Tx1East - self.Rx1East, 2) +
              np.power(Idp.Tx1North - self.Rx1North, 2) +
              np.power(Idp.Tx1Elev - self.Rx1Elev, 2)), 0.5))
        r2 = (np.power((np.power(Idp.Tx1East - self.Rx2East, 2) +
              np.power(Idp.Tx1North - self.Rx2North, 2) +
              np.power(Idp.Tx1Elev - self.Rx2Elev, 2)), 0.5))
        r3 = (np.power(np.power(Idp.Tx2East - self.Rx1East, 2) +
              np.power(Idp.Tx2North - self.Rx1North, 2), 0.5))
        r4 = (np.power(np.power(Idp.Tx2East - self.Rx2East, 2) +
              np.power(Idp.Tx2North - self.Rx2North, 2), 0.5))
        k = 2 * np.pi * (1 / ((1.0 / r1 - 1.0 / r2) - (1.0 / r3 - 1.0 / r4)))

        rho_error = (self.Vp / Idp.In *
                     np.sqrt(np.power(self.Vp_err / 100., 2) +
                             np.power(Idp.In_err / 100., 2)))

        return rho_error

    def calcRho(self, Idp):
        r1 = (np.pow((np.pow(Idp.Tx1East - self.Rx1East, 2) +
              np.pow(Idp.Tx1North - self.Rx1North, 2) +
              np.pow(Idp.Tx1Elev - self.Rx1Elev, 2)), 0.5))
        r2 = (np.pow((np.pow(Idp.Tx1East - self.Rx2East, 2) +
              np.pow(Idp.Tx1North - self.Rx2North, 2) +
              np.pow(Idp.Tx1Elev - self.Rx2Elev, 2)), 0.5))
        r3 = (np.pow((np.pow(Idp.Tx2East - self.Rx1East, 2) +
              np.pow(Idp.Tx2North - self.Rx1North, 2) +
              np.pow(Idp.Tx2Elev - self.Rx1Elev, 2)), 0.5))
        r4 = (np.pow((np.pow(Idp.Tx2East - self.Rx2East, 2) +
              np.pow(Idp.Tx2North - self.Rx2North, 2) +
              np.pow(Idp.Tx2Elev - self.Rx2Elev, 2)), 0.5))
        k = 2 * np.pi * (1 / ((1.0 / r1 - 1.0 / r2) - (1.0 / r3 - 1.0 / r4)))

        rho = k * self.Vp / Idp.In
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
        timing_string = self.headers[2].split(' ')[2]
        timing = timing_string.split(',')
        self.window_start = np.zeros(len(timing) - 1)
        self.window_end = np.zeros(len(timing) - 1)
        self.window_center = np.zeros(len(timing) - 1)
        self.window_width = np.zeros(len(timing) - 1)
        # loops each available IP window and calculate width and centre
        for i in range(len(timing) - 1):
            wintime = timing[i].split(':')
            self.window_start[i] = float(wintime[0])
            self.window_end[i] = float(wintime[1])
            self.window_center[i] = (self.window_start[i] +
                                     self.window_end[i]) / 2.0
            self.window_width[i] = (self.window_end[i] -
                                    self.window_start[i])


class Jreadtxtline:
    """
    Class specifically for reading a line of text from a dias file

    """

    def __init__(self, hdrLine, dataLine):
        # make a structure with the header inputs
        self.Vs = []
        for n in range(len(hdrLine)):
            if hdrLine[n].find("M") == 0:
                if hdrLine[n].find("MA") != 0:
                    self.Vs.append(float(dataLine[n]))
                else:
                    setattr(self, hdrLine[n], dataLine[n])
            else:
                setattr(self, hdrLine[n], dataLine[n])
