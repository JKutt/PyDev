import numpy as np
import matplotlib.pyplot as plt

###############################################
# Classes


class vsWinWeights:
    def __init__(self, num_weights):
        self.weights = np.zeros(num_weights)
        self.times = np.zeros(num_weights)

##############################################
# get average (box car equivalent)


def getWeightedWindowAve(win):
    # create time window
    timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
        2240, 2320, 2400,
        2480, 2560, 2640,
        2720, 2800, 2960,
        3120, 3280, 3440,
        3600, 3760]
    timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
        2320,2400,2480,2560,2640,2720
        ,2800,2960,3120,3280,3440
        ,3600,3760,3920]
    # find how many samples in first window
    cntr = 0
    # get time span for tail of windows
    start_tmp = timeFrom[win]
    end_tmp = timeTo[win]
    # print start_tmp, end_tmp
    for i in range(X1.size):
        if time[i] >= start_tmp and time[i] <= end_tmp:
            cntr += 1
    # create window wieghts
    indx1 = np.arange(0, 151)
    weights = 0.5 - (0.5 * np.cos((2 * np.pi * indx1) / (indx1.size - 1)))
    # plt.plot(indx1, weights, 'o-')
    # plt.show()
    # create new weights
    Wg = np.zeros(cntr)
    start_Wg = (indx1.size / 2.0 - 1.0) - (cntr / 2.0 ) + 1
    for r in range(Wg.size):
        Wg[r] = weights[int(start_Wg) + r]
    # print Wg
    # create vector storing weighted values
    Vs_window = np.zeros(cntr)
    Vs_window_ave = np.zeros(cntr)
    # get window times
    w_idx = np.zeros(cntr)
    # assign total time and time step
    count = 0
    for i in range(X1.size):
        if time[i] >= start_tmp and time[i] <= end_tmp:
            w_idx[count] = time[i]
            Vs_window[count] = Z1[i] * -1 * Wg[count]
            Vs_window_ave[count] = Z1[i] * -1
            count += 1
    #sumWin = np.sum(Vs_window)
    sumWin1 = np.sum(Vs_window_ave)
    # print Vs_window
    return sumWin1 / cntr
    #print sumWin1 / cntr
    
######################################################
# get tappered windowing


def getWeightedWindow(win, overlap, tapper):
    # create time window
    timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
        2240, 2320, 2400,
        2480, 2560, 2640,
        2720, 2800, 2960,
        3120, 3280, 3440,
        3600, 3760]
    timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
        2320,2400,2480,2560,2640,2720
        ,2800,2960,3120,3280,3440
        ,3600,3760,3920]
    # find how many samples in first window
    cntr = 0
    cnt_expand = 0
    # get time span for tail of windows
    if (win < overlap):
        end_tmp = timeTo[win + cnt_expand]
        start_tmp = timeFrom[win - cnt_expand]
        cnt_expand = cnt_expand + 1
    elif (win + overlap) >= 19:
        end_tmp = timeTo[19]
        start_tmp = timeFrom[win - overlap]
    elif (win - overlap) <= 0:
        start_tmp = timeFrom[0]
        end_tmp = timeTo[win + overlap]
    else:
        # start_tmp = timeFrom[win] - overlap * (timeTo[win - 1] - timeFrom[win - 1])
        start_tmp = timeFrom[win - overlap]
        end_tmp = timeTo[win + overlap]

        # end_tmp = timeTo[win] + overlap * (timeTo[win] - timeFrom[win])
    center = 0
    window_center = ((timeTo[win] + timeFrom[win]) / 2.0)
    for i in range(X1.size):
        if time[i] >= start_tmp and time[i] <= end_tmp:
            cntr += 1
            if time[i] > window_center and center == 0:
                center = cntr
    # print "number of sample are:"
    # create window wieghts
    indx1 = np.arange(0, tapper)
    weights = 0.5 - (0.5 * np.cos((2 * np.pi * indx1) / (indx1.size - 1)))
    # create new weights
    Wg = np.zeros(cntr)
    if win == 0:
        start_Wg = (indx1.size / 2.0 - 1.0) - (center) + 1
    elif win == 19:
        start_Wg = (indx1.size / 2.0 - 1.0) - (center) + 1
    else:
        start_Wg = (indx1.size / 2.0 - 1.0) - (center) + 1
    for r in range(Wg.size):
        Wg[r] = weights[int(start_Wg) + r]
    # print Wg
    # create vector storing weighted values
    Vs_window = np.zeros(cntr)
    Vs_window_ave = np.zeros(cntr)
    # get window times
    w_idx = np.zeros(cntr)
    # assign total time and time step
    count = 0
    for i in range(X1.size):
        if time[i] >= start_tmp and time[i] <= end_tmp:
            w_idx[count] = time[i]
            Vs_window[count] = Z1[i] * -1 * Wg[count]
            Vs_window_ave[count] = Z1[i] * -1
            count += 1
    sumWin = np.sum(Vs_window)
    # print Vs_window
    return sumWin / cntr, Wg, w_idx
    # print sumWin1 / cntr


def getWeightedWindow2(win, overlap, tapper):
    # create time window
    timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
        2240, 2320, 2400,
        2480, 2560, 2640,
        2720, 2800, 2960,
        3120, 3280, 3440,
        3600, 3760]
    timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
        2320,2400,2480,2560,2640,2720
        ,2800,2960,3120,3280,3440
        ,3600,3760,3920]
    # find how many samples in first window
    cntr = 0
    # get time span for tail of windows
    start_tmp = timeFrom[win] - overlap * (timeTo[win - 1] - timeFrom[win - 1])
    if start_tmp < timeFrom[0]:
        start_tmp = timeFrom[0]

    end_tmp = timeTo[win] + overlap * (timeTo[win] - timeFrom[win])
    if end_tmp > timeTo[19]:
        end_tmp = timeTo[19]
        # end_tmp = timeTo[win] + overlap * (timeTo[win - 1] - timeFrom[win - 1])
    # print start_tmp, end_tmp
    center = 0
    window_center = ((timeTo[win] + timeFrom[win]) / 2.0)
    for i in range(X1.size):
        if time[i] >= start_tmp and time[i] <= end_tmp:
            cntr += 1
            if time[i] > window_center and center == 0:
                center = cntr
    # print "number of sample are:"
    # print cntr
    # create window wieghts
    indx1 = np.arange(0, tapper)
    weights = 0.5 - (0.5 * np.cos((2 * np.pi * indx1) / (indx1.size - 1)))
    # create new weights
    Wg = np.zeros(cntr)
    start_Wg = (indx1.size / 2.0 - 1.0) - (center) + 1
    for r in range(Wg.size):
        Wg[r] = weights[int(start_Wg) + r]
    # print Wg
    # create vector storing weighted values
    Vs_window = np.zeros(cntr)
    Vs_window_ave = np.zeros(cntr)
    # get window times
    w_idx = np.zeros(cntr)
    # assign total time and time step
    count = 0
    for i in range(X1.size):
        if time[i] >= start_tmp and time[i] <= end_tmp:
            w_idx[count] = time[i]
            Vs_window[count] = Z1[i] * -1 * Wg[count]
            Vs_window_ave[count] = Z1[i] * -1
            count += 1
    sumWin = np.sum(Vs_window)
    # print Vs_window
    return sumWin / cntr, Wg, w_idx
#######################################################


fileName = "/Users/juan/Documents/testData/R92_J5_C0.xyz"
size_stack = 600                            # size of stacked data
sample_rate = 150                           # sample rate of data
time_step = 4.0 / 600.0
# open the file
text_file = open(fileName, "r")
lines = 0
# determin how many lines in the file
while text_file.readline():
        lines += 1
text_file.close()
# intiate numpy variables
X1 = np.zeros((lines))
Z1 = np.zeros((lines))
# initiate data count
cnt = 0
# initiate header count
cnthdr = 0
# open file again
text_file = open(fileName, "r")
# retrieve and assign the data
for line in text_file:
        spl1 = line.split()
        X1[cnt] = float(cnt)
        Z1[cnt] = float(spl1[0])
        cnt += 1
text_file.close()
amp = np.ones(20)
time = X1 * time_step * 1000

# plt.plot(time, Z1)
# plt.show()

#####################################################
# Calculate the decay
timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
        2240, 2320, 2400,
        2480, 2560, 2640,
        2720, 2800, 2960,
        3120, 3280, 3440,
        3600, 3760]
timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
        2320,2400,2480,2560,2640,2720
        ,2800,2960,3120,3280,3440
        ,3600,3760,3920]
W2 = np.zeros(20)
W2_2 = np.zeros(20)
W2_3 = np.zeros(20)
W2_4 = np.zeros(20)
W2_5 = np.zeros(20)
t_ = np.zeros(20)
Wg = []
for i in range(W2.size):
    W2[i], Wg_tmp, t_tmp = getWeightedWindow2(i, 0.45, 401)
    W2_3[i], Wg_tmp, t_tmp = getWeightedWindow2(i, 0.95, 301)
    W2_4[i], Wg_tmp, t_tmp = getWeightedWindow(i, 3, 301)
    W2_5[i], Wg_tmp, t_tmp = getWeightedWindow(i, 4, 301)
    W2_2[i] = getWeightedWindowAve(i)
    Wg_w2 = vsWinWeights(Wg_tmp.size)
    Wg_w2.weights = Wg_tmp
    Wg_w2.times = t_tmp
    Wg.append(Wg_w2)
    t_[i] = (timeFrom[i] + timeTo[i]) / 2.0
print Wg[19].weights
# plot the decay
f, plts = plt.subplots(2, 2)
plts[1, 0].plot(t_, W2, 'm-o')
plts[1, 0].plot(t_, W2_3, 'b-o')
plts[1, 0].plot(t_, W2_4, 'c-o')
plts[1, 0].plot(t_, W2_5, 'g-o')
# plts[1, 0].st_xscale("log", nonposx='clip')
# plts[1, 0].set_yscale("log", nonposy='clip')
plts[1, 0].set_ylabel('Voltage (mV)')
plts[1, 0].set_xlabel('Time (ms)')
plts[1, 0].set_title('tapper Windows')
plts[1, 1].plot(t_, W2_4, 'b-o')
plts[1, 1].set_ylabel('Voltage (mV)')
plts[1, 1].set_xlabel('Time (ms)')
plts[1, 1].set_title('Averaged')
plts[0, 1].plot(t_, W2, 'm-o')
plts[0, 1].plot(t_, W2_2, 'k-o')
plts[0, 1].plot(t_, W2_3, 'b-o')
plts[0, 1].plot(t_, W2_4, 'c-o')
plts[0, 1].set_title('Extracted Vs - 45% & 95% & MG overlap')
plts[0, 1].set_ylabel('Voltage (mV)')
plts[0, 1].legend(['win Samp #401',
                  'Ave', 'win Samp #301',
                   '3 gates tapper'])
plts[0, 1].grid()
plts[0, 0].plot(time[320:580], Z1[320:580] * -1., 'g-o')
plts[0, 0].plot(Wg[0].times, Wg[0].weights, 'm-')
plts[0, 0].plot(Wg[1].times, Wg[1].weights, 'm-')
plts[0, 0].plot(Wg[2].times, Wg[2].weights, 'm-')
plts[0, 0].plot(Wg[3].times, Wg[3].weights, 'm-')
plts[0, 0].plot(Wg[4].times, Wg[4].weights, 'm-')
plts[0, 0].plot(Wg[5].times, Wg[5].weights, 'm-')
plts[0, 0].plot(Wg[6].times, Wg[6].weights, 'm-')
plts[0, 0].plot(Wg[7].times, Wg[7].weights, 'm-')
plts[0, 0].plot(Wg[8].times, Wg[8].weights, 'm-')
plts[0, 0].plot(Wg[9].times, Wg[9].weights, 'm-')
plts[0, 0].plot(Wg[10].times, Wg[10].weights, 'm-')
plts[0, 0].plot(Wg[11].times, Wg[11].weights, 'm-')
plts[0, 0].plot(Wg[12].times, Wg[12].weights, 'm-')
plts[0, 0].plot(Wg[13].times, Wg[13].weights, 'm-')
plts[0, 0].plot(Wg[14].times, Wg[14].weights, 'm-')
plts[0, 0].plot(Wg[15].times, Wg[15].weights, 'm-')
plts[0, 0].plot(Wg[16].times, Wg[16].weights, 'm-')
plts[0, 0].plot(Wg[17].times, Wg[17].weights, 'm-')
plts[0, 0].plot(Wg[18].times, Wg[18].weights, 'm-')
plts[0, 0].plot(Wg[19].times, Wg[19].weights, 'm-')
plts[0, 0].grid()
plts[0, 0].set_ylabel('V')
plts[0, 0].set_title('Stacked Time Series')
plt.show()
##################################################
