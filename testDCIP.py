import unittest
import numpy as np
import DCIPsimpeg.DCIPtools as DCIP

class TestingReCalculationMethods(unittest.TestCase):

    # testing rho calculation
    def testRhoCalculation(self):
        # homogenous halfspace
        rho = 2200
        In = 1.00

        # create synthetic variables testing for a dipole generation and its calculations
        hdr_line = ['DIPOLE', 'RDG', 'Status', 'Rx1x', 'Rx1y', 'Rx2x', 'Rx2y', 'Tx1x', 'Tx1y', 'Tx2x', 'Tx2y',
                    'Rx1East', 'Rx1North', 'Rx1Elev', 'Rx2East', 'Rx2North', 'Rx2Elev',
                    'Tx1East', 'Tx1North', 'Tx1Elev', 'Tx2East', 'Tx2North', 'Tx2Elev',
                    'Rx1File', 'Rx2File', 'TxFile', 'k',
                    'Coupling', 'Sp', 'Vp_err', 'Vp', 'In', 'In_err', 'Rho', 'Rho_QC', 'Mx_QC',
                    'Mx_err', 'Stack', 'TimeBase', 'Vs']
        data_line = ['1', '1', 'Good', '500', '500', '700', '500', '200', '500', '100', '500',
                     '345500', '4345500', '520', '345700', '4345500', '522',
                     '345200', '4345500', '519', '345100', '4345500', '520', 'AA000.raw', 'AA000.raw',
                      'AA000.raw', '0',
                     '15.0', '0.0', '0.0', '0.0', str(In), '0.0', str(rho), 'Accept', 'Accept', '0.0',
                     '0', '2000', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0',
                     '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0']
        varInfo = DCIP.Jreadtxtline(hdr_line, data_line)
        voltage_dipole = DCIP.JvoltDipole(varInfo)
        current_dipole = DCIP.JinDipole(varInfo)
        
        # calculate expected vp of the homogenous halfspace
        gk = voltage_dipole.calcGeoFactor(current_dipole)
        vp = In * rho * (1 / gk)
        # assign this volatage for testing
        voltage_dipole.Vp = vp
        # calc the resistivity with synthetic Vp to see if it matches known
        rho_test = voltage_dipole.calcRho(current_dipole)
        # send results
        self.assertTrue(np.isclose(rho_test, rho))

    # test mx calculation
    def testMxCalculation(self):
        # homogenous halfspace
        rho = 2200
        mx = 1
        In = 1.00

        # window widths
        windows = [20,20,40,40,40,40,80,80,80,80,80,80,80,160,160,160,160,160,160,160]
        hdr_line = ['DIPOLE', 'RDG', 'Status', 'Rx1x', 'Rx1y', 'Rx2x', 'Rx2y', 'Tx1x', 'Tx1y', 'Tx2x', 'Tx2y',
                    'Rx1East', 'Rx1North', 'Rx1Elev', 'Rx2East', 'Rx2North', 'Rx2Elev',
                    'Tx1East', 'Tx1North', 'Tx1Elev', 'Tx2East', 'Tx2North', 'Tx2Elev',
                    'Rx1File', 'Rx2File', 'TxFile', 'k',
                    'Coupling', 'Sp', 'Vp_err', 'Vp', 'In', 'In_err', 'Rho', 'Rho_QC', 'Mx_QC',
                    'Mx_err', 'Stack', 'TimeBase', 'Vs']
        data_line = ['1', '1', 'Good', '500', '500', '700', '500', '200', '500', '100', '500',
                     '345500', '4345500', '520', '345700', '4345500', '522',
                     '345200', '4345500', '519', '345100', '4345500', '520', 'AA000.raw', 'AA000.raw',
                      'AA000.raw', '0',
                     '15.0', '0.0', '0.0', '0.0', str(In), '0.0', str(rho), 'Accept', 'Accept', '0.0',
                     '0', '2000', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0',
                     '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0', '0.0']
        varInfo = DCIP.Jreadtxtline(hdr_line, data_line)
        # create voltage dipole
        voltage_dipole = DCIP.JvoltDipole(varInfo)
        # create current dipole
        current_dipole = DCIP.JinDipole(varInfo)
        
        # calculate expected vp of the homogenous halfspace
        gk = voltage_dipole.calcGeoFactor(current_dipole)
        vp = In * rho * (1 / gk)
        # assign this volatage for testing
        voltage_dipole.Vp = vp
        voltage_dipole.Vs = np.ones(20) * 1e-3
        # calc the resistivity with synthetic Vp to see if it matches known
        mx_test = voltage_dipole.calcMx(windows, 0, 20)
        # send results
        self.assertTrue(np.isclose(mx_test, mx / vp))


if __name__ == '__main__':
    unittest.main()
    