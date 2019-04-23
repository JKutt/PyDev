import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
# from multiprocessing import Pool

# ============================ User area ==============================
# ver 0.1
# set variables
fname = "E:/Projects/China/Karst/IPDB/B1_NorthBlock-QC.DAT"
outname = "E:/Projects/China/Karst/IPDB/Lusha_North-QC.log"
utm_csv = "E:/Projects/China/Karst/IPDB/NorthBlock-Karst-DEM.csv"
ddn_info = "E:/Projects/debug/Gunjan/orano/mem_I.txt"

# =====================================================================
"""
   Load the DIAS database along with the handheld gps/DEM data
"""
# load the data file
patch = DCIP.loadDias(fname)
# load DEM CSV file
utms, local, level = DCIP.loadCsvDem(utm_csv)
# load DDN info
ddn = np.genfromtxt(ddn_info)


# =====================================================================
# 1) compares dem to dias database file
"""
   evaluates the handheld gps data from field

   OUTPUT: plots of local and utm locations
"""
# plots diferences in locations adjacent
DCIP.plotCsvDemLocationDifferences(local=local, utm=utms)
# plots survey grid local + utm
DCIP.plotLocalAndUtmGrids(local=local, utm=utms)


# =====================================================================
# 2) check the continuity of the files used
"""
    goes through database and creates node database for each reading
    and checks that file number id's are consecutive

    OUTPUT: log file indicating discontinuities
"""
patch.checkContinuityInNodeFiles(path=outname)


# =====================================================================
# 3) checks records in log file to records actually processed
# """
#     runs through database and compares availble records to
#     records log
# """
rec_log_fpath = "E:/Projects/China/Karst/IPDB/B2/Log/Recordings_0.txt"
patch.checkForMissingRecords(record_log=rec_log_fpath)


# =====================================================================
# 4) Plot both GPS and database locations
"""
   PLots the GPS data from the handheld over database
   locations to ensure all stations covered

   OUTPUT: plot of the overlay projection
"""
# plots all electrode locations
patch.plotGpsOverDatabaseLocations(utms, rx=True, tx=True, dipole_dipole=False)
# plots rx locations
patch.plotGpsOverDatabaseLocations(utms, rx=True, tx=False, dipole_dipole=False)
# plots tx locations
patch.plotGpsOverDatabaseLocations(utms, rx=False, tx=True, dipole_dipole=False)
# plots all electrode locations local
patch.plotGpsOverDatabaseLocations(local, rx=True, tx=True, dipole_dipole=False, local=True)
# plots rx locations local
patch.plotGpsOverDatabaseLocations(local, rx=True, tx=False, dipole_dipole=False, local=True)
# plots tx locations local
patch.plotGpsOverDatabaseLocations(local, rx=False, tx=True, dipole_dipole=False, local=True)

# # =====================================================================
# # 5) Plot Elevations
# """
#    PLots the Elevation data from the DEM/GPS of input data and database 

#    OUTPUT: plot of the elevationms
# """
# patch.plotElevations(utms, 2D=False)
# TODO !!!!!!


# =====================================================================
# 6) Current comparisons
"""
    Compares current input from field DDN to database and plots
.    percent difference and the currents from each
"""
patch.compareInjectedCurrents2fieldDDN(ddn=ddn, reject=None, line_dir='ns')
patch.plotLabeledInjections(gps=None)

# =====================================================================
# 7) plot node locations + Node ID for a defined reading
"""
    runs through each reading and plots the avaible nodes and dipole
    polarity
"""
patch.checkDipoleAppRhoPolarityPerReading(num_readings=2,
                                          # gps_locations=utms,
                                          dipole_dipole=False)


# =====================================================================
# 8) histograms of data Vitals
# """
#    runs through database and creates histogram of the data
#    vitals such as rho, Mx, Vp and prints to screen  some stats
# """
patch.getPercentageLowVoltages(threshold=1.0, reject=None)
patch.getPercentageNegativeAppMx(reject=None)
patch.getPercentageNegativeAppRho(reject=None)
patch.plotHistogramOfDataVitals(log_rho=False, log_vp=False, reject='app_rho')


# ========================================================================
# 9) compare recalculation of data to database values
patch.CompareDatabase2RecalcValues(0, 20, reject='app_mx')
