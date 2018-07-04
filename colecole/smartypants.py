import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP

fileName = "C:/Projects/springs/JK-test.DAT"
patch = DCIP.loadDias(fileName)

rho = []
mx =[]

for rdg in range(len(patch.readings)):
	for dp in range(len(patch.readings[rdg].Vdp)):
		rho.append(patch.readings[rdg].Vdp[dp].Rho)
		mx.append(patch.readings[rdg].Vdp[dp].Mx)

rho = np.asarray(np.abs(rho))
mx = np.asarray(mx)


win_length = 1001
zero_pad = (win_length - 1) / 2.
sides = zero_pad
kappa = 1.4826
eta = 3.0
reject = []
value = []
for index in range(mx.size):
	win = np.zeros(win_length)
	try:
		if (index - (win_length - 1) / 2.) <= 0:
			diff = int(np.abs((index - ((win_length - 1) / 2.0))))
			win[diff:win_length] = mx[0:(index + int(sides) + 1)]
			local_median = np.median(win)
			mean_abs_dev = np.median(np.abs(win - local_median))
			sigma = eta * mean_abs_dev
			lhs = np.abs(mx[index] - local_median)
			rhs = eta * sigma
			if lhs > rhs:
				print("found outlier")
				reject.append(index)
				value.append(mx[index])
		elif (index + (win_length - 1) / 2.) >= (mx.size - 1):
			diff = int(np.abs((index - ((win_length - 1) / 2.0))))
			win[0:(win_length - 1) - diff] = mx[(index - int(sides)) -1 :]
			local_median = np.median(win)
			mean_abs_dev = np.median(np.abs(win - local_median))
			sigma = eta * mean_abs_dev
			lhs = np.abs(mx[index] - local_median)
			rhs = eta * sigma
			if lhs > rhs:
				print("found outlier") 
				reject.append(index)
				value.append(mx[index])
		else:
			sides = int(sides)
			win = mx[(index - sides):(index + sides)]
			local_median = np.median(win)
			mean_abs_dev = np.median(np.abs(win - local_median))
			sigma = eta * mean_abs_dev
			lhs = np.abs(mx[index] - local_median)
			rhs = eta * sigma
			if lhs > rhs:
				print("found outlier") 
				reject.append(index)
				value.append(mx[index])
	except:
		pass

plt.plot(mx, '.')
plt.plot(reject, value, 'or')
plt.show()
