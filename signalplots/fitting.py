import matplotlib.pyplot as plt
import numpy as np
import Jtools as jt
import scipy.interpolate as inter

fileName = "C:\\Projects\\synctest2.xyz"

lines = 0;
splarray = " "

# open the file
text_file = open(fileName,"r")

#determin how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()
# intiate numpy variables
X = np.zeros((lines-1))
Z = np.zeros((lines-1))

# initiate data count
cnt = 0
# initiate header count
cnthdr = 0
# open file again
text_file = open(fileName,"r")
# setup counter for offtimes
offcnt = 0
cntoff = 0
# retrieve and assign the data
for line in text_file:
    if cnthdr==0:
        #print(line)
        cnthdr += 1
    else:
        spl = line.split()
        X[cnt] = float(spl[0])
        Z[cnt] = float(spl[1])
        cnt += 1
        if offcnt > 300:
            if offcnt < 550 and offcnt > 350:
                cntoff +=1

            if offcnt == 600:
                offcnt = 0
            else:
                offcnt +=1
        else:
            offcnt +=1
        
text_file.close()

# create offtime data variable
xoff = np.zeros(cntoff)
yoff = np.zeros(cntoff)
offcnt = 0
cntrr = 0
for i in range(X.size):
    if offcnt > 300:
        if offcnt < 550 and offcnt > 350:
            xoff[cntrr] = X[i]
            yoff[cntrr] = Z[i]
            cntrr +=1

        if offcnt == 600:
            offcnt = 0
        else:
            offcnt +=1
    else:
        offcnt +=1

# try poly fit
z = np.polyfit(xoff,yoff,3)
x_new = np.linspace(X[0], X[-1], X.size)
f = np.poly1d(z)
print(z)
# calculate new x's and y's
#x_new = np.linspace(X[0], X[-1], X.size)
y_new = f(x_new)
print(x_new.size)
#stackcheda = jt.stack(Z,0.125,150)
#s1 = inter.InterpolatedUnivariateSpline (X, Z)
s1 = inter.UnivariateSpline(xoff, yoff)
xs = np.linspace(xoff[0], xoff[xoff.size-1], 20000)
#stack = np.zeros([600])
#print(stack.size)
#cntr = 0
#for i in range(600):
#        stack[cntr] = Z[i]
#        cntr+= 1
##print(s1)
#plt.plot(stackcheda,'green')
# plot the data
plt.plot(X,Z-y_new,'g')
##plt.ylabel('signal')
plt.hold(True)
#plt.plot(x_new,y_new,'r')
#plt.plot(xoff,yoff,'b')
#plt.plot(xs,s1(xs),'r')
##cntr = 0
##for i in range(Z.size-600):
##        stack[cntr] = abs(Z[i])
##        cntr+= 1
##        if cntr == 600:
##                plt.plot(stack)
##                cntr = 0

plt.grid(True)
plt.show()

print("Plotting Complete")
