import matplotlib.pyplot as plt
import numpy as np
import Jtools as jt
import scipy.interpolate as inter

fileName = "C:\\Users\\johnk\\Projects\\hampel.xyz"
fileName2 = "C:\\Users\\johnk\\Projects\\hampel-applied.xyz"

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
    Z[cnt] = float(spl[1])
    cnt += 1

text_file.close()

# ======================================
lines = 0
splarray = " "

# open the file
text_file = open(fileName2, "r")

#determin how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
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
    Z2[cnt] = float(spl[1])
    cnt += 1

text_file.close()

# for i in range(X.size - 1):
#     if X[i] > X[i + 1]:
#         print("found messup at: ", i)

# take a look at the sample widths from raw
# X = (X - X[0]) * 86400.
# X2 = (X2 - X2[0]) * 86400.
# sample_width = X2[:(X2.size - 1)] - X2[1:]
# sample_width2 = X[:(X.size - 1)] - X[1:]
# print(sample_width)
# plot the data
# plt.plot(np.abs(sample_width), '.r')
# plt.plot(np.abs(sample_width2), '.g')
plt.title("Hampel Filter ( 1 sigma reject & 9 point window )")
plt.xlabel("sample")
plt.ylabel("Voltage (mV)")
# plt.plot(X, (Z - Z2) - 300., 'go-')
plt.plot(X, Z, 'ko-')
plt.plot(X2, Z2, 'ro-')
# plt.plot(X2[9299], Z2[9299], 'dm')
# plt.plot(X[9299], Z[9299], 'dc')
# for i in range(Z.size-600):
#        stack[cntr] = abs(Z[i])
#        cntr+= 1
#        if cntr == 600:
#                plt.plot(stack)
#                cntr = 0

plt.grid(True)
plt.show()

print("Plotting Complete")
