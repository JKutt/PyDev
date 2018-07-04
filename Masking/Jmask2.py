import matplotlib.pyplot as plt
import numpy as np

fileName = "C:/Users/johnk/devProjects/Python/Masking/L0_MxA.xyz"
# outName = "C:/Users/johnk/devProjects/Python/Masking/L0_MxA.csv"
outName = "C:/Users/johnk/devProjects/Python/Masking/L0_MxA_mask.ply"

lines = 0
splarray = " "

# open the file
text_file = open(fileName, "r")

# determin how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()
# intiate numpy variables
X = np.zeros((lines - 1))
Z = np.zeros((lines - 1))

# initiate data count
cnt = 0
# initiate header count
cnthdr = 0
# open file again
text_file = open(fileName, "r")
# retrieve and assign the data
for line in text_file:
    if cnthdr == 0:
        # print(line)
        cnthdr += 1
    else:
        spl = line.split()
        X[cnt] = float(spl[3])
        Z[cnt] = float(spl[4])
        cnt += 1

text_file.close()

# click point
fig = plt.figure()
ax = fig.add_subplot(111)
# ax.set_xlim([0, 10])
# ax.set_ylim([0, 10])


def onclick(event):
    # write to file
    out_file = open(outName, "a")
    out_file.write("%s, %s\n" % (event.xdata, event.ydata))
    out_file.close()
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    plt.plot(event.xdata, event.ydata, 'o')
    fig.canvas.draw()


cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.plot(X, Z, 'g^')
plt.show()

print("Masking Complete")
