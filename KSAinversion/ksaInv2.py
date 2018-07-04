from SimPEG import DC, IP
from SimPEG import Maps, Utils
from SimPEG import Mesh
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from time import clock
from pylab import hist
import DCIPtools as tools
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver

filename = "/Users/juan/Documents/testData/KSA/3D_QC_Preliminary.DAT"
filename1 = "/Users/juan/Documents/testData/KSA/ksa3D-DC.con"
filenameObs = "/Users/juan/Documents/testData/KSA/ksa3D-IP.obs"
fileName2 = "/Users/juan/Documents/testData/Mesh.msh"
mesh = Mesh.TensorMesh._readUBC_3DMesh(fileName2)
patch = tools.loadDias(filename)
print(len(patch.readings))
survey_dc, tx = patch.createDcSurvey("DC")
survey_dc.survey_type = "dipole-dipole"
survey_dc.data_dc_type = 'volt'
survey_dc.getABMN_locations()

# uniq = Utils.uniqueRows(np.vstack((survey_dc.a_locations,
#                                    survey_dc.b_locations,
#                                    survey_dc.m_locations,
#                                    survey_dc.n_locations)))
# electrode_locations = uniq[0]                            # assign
# actinds = Utils.surface2ind_topo(mesh,
#                                  electrode_locations,
#                                  method='cubic')      # active indicies
# survey_dc.drapeTopo(mesh, actinds)

survey_ip1, tx = patch.createDcSurvey("IP")
plt.plot(survey_ip1.dobs, '.')
plt.show()
