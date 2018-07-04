import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP


def getTime():
    timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
                2240, 2320, 2400,
                2480, 2560, 2640,
                2720, 2800, 2960,
                3120, 3280, 3440,
                3600, 3760]
    timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
              2320, 2400, 2480, 2560, 2640, 2720,
              2800, 2960, 3120, 3280, 3440,
              3600, 3760, 3920]
    return timeFrom, timeTo


fileName = "/home/juan/Documents/DevProjects/testData/L5150-final.DAT"
patch = DCIP.loadDias(fileName)
timeFrom, timeTo = getTime()
timeTo = np.asarray(timeTo)
timeFrom = np.asarray(timeFrom)
timeCenter = (timeTo + timeFrom) / 2.
for go in range(len(patch.readings[0].Vdp)):
	decay = (patch.readings[0].Vdp[go].Vs /
	         patch.readings[0].Vdp[go].Vp)
	w_ = np.ones(timeCenter.size)             # cole-cole weights
	w_[:3] = 0.3
	c = 0.65
	tau = 0.35
	r = 6.0
	stored_error = 0
	min_iter = 3
	num_windows = timeTo.size

	# for iter in range(12):
	#     c, tau, M, error, vs = DCIP.getColeCole(decay,
	#                                             c,
	#                                             tau,
	#                                             r,
	#                                             timeCenter,
	#                                             w_)
	#     delta_error = abs(stored_error - error) / ((stored_error + error) / 2.)
	#     print("iter: %i | c: %f | tau: %e | M: %f | error: %f | delta: %f" %
	#           (iter, c, tau, M, error, delta_error))
	#     stored_error = error
	#     # r = r / 2.
	#     if delta_error > 0.002 or iter < min_iter:
	#         r = r / 2.
	#     elif delta_error < 0.002:
	#         print("convergence accomplished! DONE")
	#         break

	# print(c, tau, M, error)
	# percent_diff = (
	#     np.sum((np.abs((decay - vs)) /
	#            ((vs + decay) / 2.)) * w_)) / num_windows
	# global_mean = np.abs((decay - vs) /
	#            ((vs + decay) / 2.) * w_)
	# print(percent_diff)
	# print("mean diff: ", np.abs((decay - vs)) /
	#            ((vs + decay) / 2.) * w_)
	# print("mean diff Global", np.mean(global_mean))
	# print("median diff Global", np.median(global_mean))
	# print("median diff Global", np.std(global_mean))
	# print("diff: ", np.abs((decay - vs)))
	# print(decay)

	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

	ax1.plot(timeCenter, decay, 'o-')
	# ax1.plot(timeCenter, vs, 'o-r')
	# ax2.hist(global_mean, 10)
	# ax3.plot(timeCenter, global_mean, 'o-')

	# for i in range(1):
	#     for j in range(len(patch.readings[i].Vdp)):
	#         decay = patch.readings[i].Vdp[j].Vs
	#         plt.plot(timeCenter, decay, 'o-')
	# plt.show()


	# rho = []
	# mx =[]

	# for rdg in range(len(patch.readings)):
	# 	for dp in range(len(patch.readings[rdg].Vdp)):
	# 		rho.append(patch.readings[rdg].Vdp[dp].Rho)
	# 		mx.append(patch.readings[rdg].Vdp[dp].Mx)

	# rho = np.asarray(np.abs(rho))
	# mx = np.asarray(mx)
	deriv = np.diff(decay)
	ax4.plot(timeCenter[1:], deriv, '-ob')
	ax1.grid(True)
	ax1.plot(timeCenter[1:], deriv, '-ob')

	win_length = 7
	zero_pad = (win_length - 1) / 2.
	sides = zero_pad
	kappa = 1.4826
	eta = 1.5
	reject = []
	value = []
	positive_deriv = 0
	for index in range(deriv.size):
		win = np.zeros(win_length)
		if deriv[index] >= 0:
			positive_deriv = positive_deriv + 1
			reject.append(timeCenter[index + 1])
			value.append(deriv[index])
		# try:
		# 	# if global_mean[index] > 0.3:
		# 	# 	print("found Outlier") 
		# 	# 	reject.append(timeCenter[index])
		# 	# 	value.append(global_mean[index])
		# 	if (index - (win_length - 1) / 2.) <= 0:
		# 		diff = int(np.abs((index - ((win_length - 1) / 2.0))))
		# 		win[diff:win_length] = deriv[0:(index + int(sides) + 1)]
		# 		local_median = np.median(win)
		# 		mean_abs_dev = np.median(np.abs(win - local_median))
		# 		sigma = eta * mean_abs_dev
		# 		lhs = np.abs(deriv[index] - local_median)
		# 		rhs = eta * sigma
		# 		# print("lhs: ", lhs, " rhs: ", rhs)
		# 		# print("mad: ", global_mean, " lm: ", local_median)
		# 		if lhs > rhs:
		# 			print("found outlier")
		# 			reject.append(timeCenter[index + 1])
		# 			value.append(deriv[index])
		# 	elif (index + (win_length - 1) / 2.) >= (deriv.size - 1):
		# 		diff = int(np.abs((index - ((win_length - 1) / 2.0))))
		# 		win[0:(win_length - 1) - diff] = deriv[(index - int(sides)) -1 :]
		# 		local_median = np.median(win)
		# 		mean_abs_dev = np.median(np.abs(win - local_median))
		# 		sigma = eta * mean_abs_dev
		# 		lhs = np.abs(deriv[index] - local_median)
		# 		rhs = eta * sigma
		# 		# print("lhs: ", lhs, " rhs: ", rhs)
		# 		# print("mad: ", global_mean, " lm: ", local_median)
		# 		if lhs > rhs:
		# 			print("found outlier") 
		# 			reject.append(timeCenter[index + 1])
		# 			value.append(deriv[index])
		# 	else:
		# 		sides = int(sides)
		# 		win = deriv[(index - sides):(index + sides)]
		# 		local_median = np.median(win)
		# 		mean_abs_dev = np.median(np.abs(win - local_median))
		# 		sigma = eta * mean_abs_dev
		# 		lhs = np.abs(deriv[index] - local_median)
		# 		rhs = eta * sigma
		# 		# print("lhs: ", lhs, " rhs: ", rhs)
		# 		# print("mad: ", global_mean, " lm: ", local_median)
		# 		if lhs > rhs:
		# 			print("found outlier") 
		# 			reject.append(timeCenter[index + 1])
		# 			value.append(deriv[index])
		# except:
		# 	pass
	if len(value) > (0.2 * deriv.size):
		print("rejected")
	# if positive_deriv > 2:
	# 	print("rejected bc positive")
	# ax4.plot(timeCenter, global_mean, '.')
	ax4.plot(reject, value, 'or')

	plt.show()
