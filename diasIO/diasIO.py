import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP

# define the file required for import
fileName = "C:\\Projects\\FieldSchool\\FS2017\\FS2017-TDB.DAT"
stack_dir = "C:\\Projects\\FieldSchool\\FS2017\\Stack\\"
in_stack_dir = "C:\\Projects\\FieldSchool\\FS2017\\InStack\\"
plotIt = True
patch = DCIP.loadDias(fileName)   # Create a patch from data file

rdg = 10                           # Tx Source to plot
dp = 802                            # dipole to plot
try:
	dp_stack = patch.readings[rdg].Vdp[dp].getDipoleStack(stack_dir) # get dipole Stack
	vs_decay = patch.readings[rdg].Vdp[dp].Vs                        # this returns the windowed data
	times = patch.window_center                                      # this gets the window centers
except:
	print("corrupt Rx file")
	plotIt = False
try:
	in_stack = patch.readings[rdg].Idp.getTxStack(in_stack_dir)     # get In Stack
except:
	print("corrupt Tx file")
	plotIt = False

# Plotting time to see the requested data
if plotIt:
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	ax1.plot(patch.window_center,
	         patch.readings[rdg].Vdp[0].Vs, '-o')
	ax1.set_ylabel("Vs (mV)")
	ax1.set_xlabel("time (ms)")
	ax1.set_title("Vs Decay")
	ax2.plot(dp_stack, '-o')
	ax2.set_title("Vp stack")
	ax2.set_ylabel("Voltage (mV)")
	ax2.set_xlabel("sample #")
	ax2.grid()
	ax3.plot(in_stack, '-o')
	ax3.set_title("In stack")
	ax3.set_ylabel("Current (mA)")
	ax3.set_xlabel("sample #")
	ax3.grid()
	plt.show()
