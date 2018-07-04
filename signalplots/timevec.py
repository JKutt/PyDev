import matplotlib.pyplot as plt
import numpy as np
import Jtools as jt
import scipy.interpolate as inter

fileName = "C:\\Users\\johnk\\Projects\\test\\DataShift\\test-pre.xyz"
fileName2 = "C:\\Users\\johnk\\Projects\\test\\DataShift\\test-pre3.xyz"

lines = 0
splarray = " "

# open the file
text_file = open(fileName, "r")

# determin how many lines in the file
while text_file.readline():
        lines += 1
# print(lines)
text_file.close()
# intiate numpy variables
X = np.zeros((lines))
Z = np.zeros((lines))

# initiate data count
cnt = 0
# initiate header count
cnthdr = 0
# open file again
text_file = open(fileName, "r")
# retrieve and assign the data
for line in text_file:
    spl = line.split()
    X[cnt] = float(spl[0])
    # Z[cnt] = float(spl[1])
    cnt += 1

text_file.close()

# ======================================
lines = 0
splarray = " "

# open the file
text_file = open(fileName2, "r")

# determin how many lines in the file
while text_file.readline():
        lines += 1
# print(lines)
text_file.close()
# intiate numpy variables
X2 = np.zeros((lines))
Z2 = np.zeros((lines))

# initiate data count
cnt = 0
# initiate header count
cnthdr = 0
# open file again
text_file = open(fileName2, "r")
# retrieve and assign the data
for line in text_file:
    spl = line.split()
    X2[cnt] = float(spl[0])
    # Z2[cnt] = float(spl[1])
    cnt += 1

text_file.close()

#  =======================================
#  figure the problem
index = []
index2 = []
for i in range(X.size):
    if X[i] > 0:
        index.append(i)

for i in range(X2.size):
    if X2[i] > 0:
        index2.append(X2[i])
index2 = np.asarray(index2)
index2 = (index2 - index2[0]) * 86400.
X_ = index2 > 0.
print(index2[X_])
plt.plot(index2[X_], '.')
plt.show()
